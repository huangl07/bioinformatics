#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$con,$vcf,$dsh,$out,$ref,$anno,$dict);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"con:s"=>\$con,
	"ref:s"=>\$ref,
	"dict:s"=>\$dict,
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
open SH,">$dsh/10.annovar1.sh";
open ANNOLIST,">$out/anno.list";
open In,$vcf;
open VCF,">$out/vcf.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$vcfs)=split(/\s+/,$_);
	if (!-f $vcfs) {
		die "check $vcfs!";
	}
	print SH "java -jar /mnt/ilustre/users/dna/.env//bin//snpEff.jar -v ref -csvStats $out/$id.anno.csv -c $con $vcfs > $out/$id.anno.primary.vcf && ";
	print SH "vcftools --vcf $out/$id.anno.primary.vcf --recode-INFO ANN --recode-INFO LOF --recode-INFO NMD --recode --out $out/$id \n";
	print ANNOLIST "$id\t$out/$id.recode.vcf\t$out/$id.anno.csv\n";
	print VCF "$out/$id.recode.vcf";
}
close In;
close SH;
close VCF,
close ANNOLIST;
open Out,">$out/CombinesVCF.json\n";
print Out "{\n";
print Out "\"Combines.combines.VCFlist\": \"$out/vcf.list\",\n";
print Out "\"Combines.combines.Refdict\": \"$dict\",\n";
print Out "\"Combines.combines.Refindex\": \"$ref.fai\",\n";
print Out "\"Combines.combines.workdir\": \"$out\",\n";
print Out "\"Combines.combines.RefFasta\": \"$ref\"\n";
print Out "}\n";
close Out;

open SH,">$dsh/10.annovar2.sh";
print SH "cd $out/ && java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/CombinesVCF.wdl -i $out/CombinesVCF.json &&  ";
print SH "perl $Bin/bin/anno-count.pl -snp $out/snp.gene.txts -indel $out/indel.genes.txt -anno $anno -out $out/pop && ";
print SH "Rscript --input $out/pop.kegg.stat --output $out/pop.kegg --top 1&& ";
print SH "Rscript --input $out/pop.go.stat --output $out/pop.go --top 1 && ";
print SH "Rscript --input $out/pop.eggnog.stat --output $out/pop.eggnog && ";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=50G -CPU 1 --maxjob $proc $dsh/10.annovar1.sh";
#`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=50G -CPU 1 --maxjob $proc $dsh/10.annovar2.sh";
#`$job`;

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
