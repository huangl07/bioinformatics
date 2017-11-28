#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$vcf,$dOut,$dShell,$ref,$dict);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"proc:s"=>\$proc,
	"vcf:s"=>\$vcf,
	"dict:s"=>\$dict,
	"ref:s"=>\$ref,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($vcf and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$vcf=ABSOLUTE_DIR($vcf);
$ref=ABSOLUTE_DIR($ref);
open JSON,">$dOut/filter.json";
print JSON "{\n";
print JSON "\"SelectVariant.rawvcf\": \"$vcf\",\n";
print JSON "\"SelectVariant.workdir\": \"$dOut\",\n";
print JSON "\"SelectVariant.refIndex\": \"$ref.fai\",\n";
print JSON "\"SelectVariant.refFasta\": \"$ref\",\n";
print JSON "\"SelectVariant.refDict\": \"$dict\"\n";
print JSON "}\n";
close JSON;



open SH,">$dShell/09.vcf-filter.sh";
print SH "cd $dOut/ && java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/SNPFilter.wdl -i $dOut/filter.json \n ";
print SH "cd $dOut/ && java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/INDELFilter.wdl -i $dOut/filter.json \n ";
close SH;
open List,">$dOut/vcf.list";
print List "snp\t$dOut/pop.snp.filter.recode.vcf\n";
print List "indel\t$dOut/pop.indel.filter.recode.vcf\n";
close List;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=20G --CPU 8  --maxjob $proc $dShell/09.vcf-filter.sh";
print $job;
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
	fq thanslate to ref format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input genome vcf file
  -ref	<file>	input reference file
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	work sh
  -h         Help

USAGE
        print $usage;
        exit;
}
