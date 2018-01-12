#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$pop,$noref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"pop:s"=>\$pop,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
	"noref"=>\$noref,
			) or &USAGE;
&USAGE unless ($vcf and $out and $pop);
$vcf=ABSOLUTE_DIR($vcf);
$pop=ABSOLUTE_DIR($pop);

mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
$mis||=0.3;
$maf||=0.05;
$dep||=2;
$mis=1-$mis;
open SH,">$dsh/01.filtered.sh";
print SH "vcftools --vcf $vcf  --out $out/pop --max-missing $mis --maf $maf --minDP $dep --recode \n ";
close SH;
if ($noref) {
	open SH,">$dsh/02.calculate.sh";
	open SH2,">$dsh/03.draw-select.sh";
	open In,$pop;
	my %group;
	my %filehand;
	my %pophand;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$gid)=split(/\s+/,$_);
		if (!exists $filehand{$gid}) {
			open $filehand{$gid},">$out/$gid.list";
		}
		if (!exists $pophand{$gid}) {
			open $pophand{$gid},">$out/$gid.pophand";
		}
		print {$filehand{$gid}} "$id\n";
		print {$pophand{$gid}} "$id\t$gid\n";
	}
	close In;
	my @groid=sort keys %filehand;
	for (my $i=0;$i<@groid;$i++) {
		print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i]  --window-pi 2000000 --window-pi-step 10000 \n";
		print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i] --TajimaD 10000 \n";
		print SH "cd $out && RAiSD -n $groid[$i] -I $out/pop.recode.vcf -f -S $out/$groid[$i].list -T 50000 -k 0.005 && ";
		print SH "perl $Bin/bin/sweep_result.pl -i $out/RAiSD_Report.$groid[$i] -o $out/$groid[$i].raisd.out\n";
		print SH2 "Rscript $Bin/bin/pi-tajima.R --tajima $out/$groid[$i].Tajima.D --pi $out/$groid[$i].windowed.pi --out $out/$groid[$i]\n";
		for (my $j=$i+1;$j<@groid;$j++) {
			open SPID,">$out/$groid[$i]-$groid[$j].pid";
			print SPID "PARSER_FORMAT=VCF\n";
			print SPID "VCF_PARSER_QUAL_QUESTION=\n";
			print SPID "VCF_PARSER_POP_FILE_QUESTION=$out/$groid[$i]-$groid[$j].pophand\n";
			print SPID "VCF_PARSER_PLOIDY_QUESTION=DIPLOID\n";
			print SPID "VCF_PARSER_POP_QUESTION=true\n";
			print SPID "VCF_PARSER_GTQUAL_QUESTION=\n";
			print SPID "VCF_PARSER_MONOMORPHIC_QUESTION=false\n";
			print SPID "VCF_PARSER_IND_QUESTION=\n";
			print SPID "VCF_PARSER_REGION_QUESTION=\n";
			print SPID "VCF_PARSER_READ_QUESTION=\n";
			print SPID "VCF_PARSER_PL_QUESTION=false\n";
			print SPID "VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false\n";
			print SPID "WRITER_FORMAT=GESTE_BAYE_SCAN\n";
			print SPID "GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP\n";
			close SPID; 
			print SH "cat $out/$groid[$i].pophand $out/$groid[$j].pophand > $out/$groid[$i]-$groid[$j].pophand && ";
			print SH "java -jar /mnt/ilustre/users/dna/.env//bin//PGDSpider2-cli.jar -inputformat vcf -outputformat GESTE_BAYE_SCAN -spid $out/$groid[$i]-$groid[$j].pid -inputfile $out/pop.recode.vcf -outputfile $out/$groid[$i]-$groid[$j].bayscan &&  ";
			print SH "bayescan_2.1 $out/$groid[$i]-$groid[$j].bayscan -snp -od $out -o pop  -threads 8  -pr_odds 10 -out_pilot  -out_freq && ";
			print SH "Rscript $Bin/bin/bayes.R --input $out/pop_fst.txt --outfile $out/pop.bayes\n";
			print SH "vcftools --vcf $out/pop.recode.vcf --weir-fst-pop $out/$groid[$i].list --weir-fst-pop $out/$groid[$j].list --out $out/$groid[$i]-$groid[$j] --fst-window-size 2000000 --fst-window-step 10000 \n";
			print SH2 "Rscript $Bin/bin/fst-pi.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] \n";
			print SH2 "Rscript $Bin/bin/manhattan.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] --tajima1 $out/$groid[$i].Tajima.D --tajima2 $out/$groid[$j].Tajima.D \n";
		}
	}
	close SH;
	close SH2;
}else{
	open SH,">$dsh/02.calculate.sh";
	open SH2,">$dsh/03.draw-select.sh";
	open In,$pop;
	my %group;
	my %filehand;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$gid)=split(/\s+/,$_);
		if (!exists $filehand{$gid}) {
			open $filehand{$gid},">$out/$gid.list";
		}
		print {$filehand{$gid}} "$id\n";
	}
	close In;
	my @groid=sort keys %filehand;
	for (my $i=0;$i<@groid;$i++) {
		print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i]  --window-pi 2000000 --window-pi-step 10000 \n";
		print SH "vcftools --vcf $out/pop.recode.vcf --keep $out/$groid[$i].list --out $out/$groid[$i] --TajimaD 10000 \n";
		print SH "cd $out && RAiSD -n $groid[$i] -I $out/pop.recode.vcf -f -S $out/$groid[$i].list -T 50000 -k 0.005 && ";
		print SH "perl $Bin/bin/sweep_result.pl -i $out/RAiSD_Report.$groid[$i] -o $out/$groid[$i].raisd.out\n";
		print SH "cd $out && OmegaPlus-F -input $out/pop.recode.vcf -name $groid[$i] -sampleList $out/$groid[$i].list -grid 100000 -minwin 10000 -maxwin 2000000 -seed 1000000 -threads 8 && ";
		print SH "perl $Bin/bin/sweep_result.pl -i $out/OmegaPlus_Report.$groid[$i] -o $out/$groid[$i].omega.out\n";
		print SH "cd $out && SweeD-MPFR-P -input $out/pop.recode.vcf -name $groid[$i] -sampleList $out/$groid[$i] -grid 100000  -threads 8 \n";
		print SH "perl $Bin/bin/sweep_result.pl -i $out/SweeD_Report.$groid[$i] -o $out/$groid[$i].sweed.out\n";
		print SH2 "Rscript $Bin/bin/pi-tajima.R --tajima $out/$groid[$i].Tajima.D --pi $out/$groid[$i].windowed.pi --out $out/$groid[$i]\n";
		print SH2 "Rscript $Bin/bin/sweep-manhattan.R --input $out/$groid[$i].raisd.out --out $out/$groid[$i].raisd\n";
		print SH2 "Rscript $Bin/bin/sweep-manhattan.R --input $out/$groid[$i].sweed.out --out $out/$groid[$i].sweed\n";
		print SH2 "Rscript $Bin/bin/sweep-manhattan.R --input $out/$groid[$i].omega.out --out $out/$groid[$i].omega\n";
		for (my $j=$i+1;$j<@groid;$j++) {
			print SH "vcftools --vcf $out/pop.recode.vcf --weir-fst-pop $out/$groid[$i].list --weir-fst-pop $out/$groid[$j].list --out $out/$groid[$i]-$groid[$j] --fst-window-size 2000000 --fst-window-step 10000 \n";
			print SH2 "Rscript $Bin/bin/fst-pi.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] \n";
			print SH2 "Rscript $Bin/bin/manhattan.R --fst $out/$groid[$i]-$groid[$j].windowed.weir.fst --pi1 $out/$groid[$i].windowed.pi --pi2 $out/$groid[$j].windowed.pi --out $out/$groid[$i]-$groid[$j] --tajima1 $out/$groid[$i].Tajima.D --tajima2 $out/$groid[$j].Tajima.D \n";
		}
	}
	close SH;
	close SH2;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:

  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -pop	<str>	group list
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
