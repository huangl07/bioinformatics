#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$gvcflist,$dOut,$dShell,$ref,$dict);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gvcf:s"=>\$gvcflist,
	"ref:s"=>\$ref,
	"dict:s"=>\$dict,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell
			) or &USAGE;
&USAGE unless ($gvcflist and $dOut and $dShell and $dict);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$gvcflist=ABSOLUTE_DIR($gvcflist);
$ref=ABSOLUTE_DIR($ref);
$dict=ABSOLUTE_DIR($dict);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open In,$gvcflist;
open List,">$dOut/total.gvcf.list";
open SH,">$dShell/08.gvcf-typing1.sh";
my $vcf;
my $number=0;
my %handfile;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$gvcf)=split(/\s+/,$_);
	if (!-f $gvcf) {
		die "check $gvcf!";
	}
	$number++;
	my $handfile=$number % 20;
	if (!exists $handfile{$handfile}) {
		open $handfile{$handfile},">$handfile.gvcf.list";
		open Out,">$dOut/$handfile.json\n";
		print Out "{\n";
		print Out "\"CombinesGVCF.CombineVCF.inputVCFs\": \"$dOut/$handfile.gvcf.list\",\n";
		print Out "\"CombinesGVCF.CombineVCF.Refdict\": \"$dict\",\n";
		print Out "\"CombinesGVCF.CombineVCF.Refindex\": \"$ref.fai\",\n";
		print Out "\"CombinesGVCF.CombineVCF.workdir\": \"$dOut\",\n";
		print Out "\"CombinesGVCF.CombineVCF.RefFasta\": \"$ref\",\n";
		print Out "\"CombinesGVCF.CombineVCF.Output\": \"$handfile\"\n";
		print Out "}\n";
		close Out;
		print List "$dOut/$handfile.gvcf.list","\n";
		print SH "cd $dOut/ && java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/CombineGVCF.wdl -i $dOut/$handfile.json \n";
	}
	print {$handfile{$handfile}} "$gvcf\n";
}
close List;
close SH;
close Out;
open Out,">$dOut/Gvcftyping.json\n";
print Out "{\n";
print Out "\"Gvcftyping.gvcftyping.inputVCFs\": \"$dOut/total.gvcf.list\",\n";
print Out "\"Gvcftyping.gvcftyping.Refdict\": \"$dict\",\n";
print Out "\"Gvcftyping.gvcftyping.Refindex\": \"$ref.fai\",\n";
print Out "\"Gvcftyping.gvcftyping.workdir\": \"$dOut\",\n";
print Out "\"Gvcftyping.gvcftyping.RefFasta\": \"$ref\"\n";
print Out "}\n";
close Out;
open SH,">$dShell/08.gvcf-typing2.sh";
print SH "cd $dOut/ && java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/GVCFtyping.wdl -i $dOut/Gvcftyping.json && ";
print SH "bcftools annotate --set-id +\'\%CHROM\\_\%POS\' $dOut/pop.noid.vcf -o $dOut/pop.variant.vcf\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=35G --CPU 8 --maxjob $proc $dShell/08.gvcf-typing.sh";
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
  -gvcf	<file>	input bamlist file
  -ref	<file>	input reference file
  -dict	<file>	input reference dict
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
