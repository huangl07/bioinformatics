#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$con,$vcf,$dsh,$out,$ref,$anno);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"con:s"=>\$con,
	"ref:s"=>\$ref,
	"anno:s"=>\$anno,
	"vcf:s"=>\$vcf,
	"dsh:s"=>\$dsh,
	"out:s"=>\$out,
	"proc:s"=>\$proc
			) or &USAGE;
&USAGE unless ($vcf and $con and $out);
$dsh||=$out;
$proc||=20;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$con=ABSOLUTE_DIR($con);
$vcf=ABSOLUTE_DIR($vcf);
$dsh=ABSOLUTE_DIR($dsh);
my $data="$out/data/ref";
open SH,">$dsh/09.annovar1.sh";
open ANNOLIST,">$out/anno.list";
open In,$vcf;
my $VCF="";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$vcfs)=split(/\s+/,$_);
	if (!-f $vcfs) {
		die "check $vcfs!";
	}
	print SH "java -Xmx50G -jar /mnt/ilustre/users/dna/.env//bin//snpEff.jar -v ref -csvStats $out/$id.anno.csv -c $con $vcfs > $out/$id.anno.primary.vcf && ";
	print ANNOLIST "$id\t$out/$id.anno.primary.vcf\t$out/$id.anno.csv\n";
	$VCF.=" -V $out/$id.anno.primary.vcf "
}
close In;
close SH;
close ANNOLIST;
open SH,">$dsh/10.annovar2.sh";
print SH "java -Djava.io.tmpdir=$out/tmp/ -Xmx50G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T CombineVariants -R $ref $VCF -o $out/pop.final.vcf.gz --genotypemergeoption UNSORTED -log $out/pop.merge.log && ";
print SH "perl $Bin/bin/anno-count.pl -snp $out/snp.anno.genes.txt -indel $out/indel.anno.genes.txt -anno $anno -out $out/pop && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/pop.kegg.stat --output $out/pop.kegg --top 1&& ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/pop.go.stat --output $out/pop.go --top 1 && ";
print SH "Rscript $Bin/bin/eff-enrich.R --input $out/pop.eggnog.stat --output $out/pop.eggnog --eggnog && ";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=50G -CPU 1 --maxjob $proc $dsh/09.annovar1.sh";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=50G -CPU 1 --maxjob $proc $dsh/09.annovar2.sh";
`$job`;

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
  -con	<file>	input snp eff config file
  -proc <num>   number of process for qsub,default 20
  -ref	<file>	reference fa file
  -vcf	<file>	input vcf list
  -out	<dir>	output dir
  -dsh	<dir>	output work_sh dir
  -h         Help

USAGE
        print $usage;
        exit;
}
