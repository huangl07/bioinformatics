#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dOut,$proc,$dShell);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell
			) or &USAGE;
&USAGE unless ($bamlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$bamlist=ABSOLUTE_DIR($bamlist);
$proc||=20;
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open SH,">$dShell/05.bam-mkdup.sh";
open In,$bamlist;
open Out,">$dOut/bam.mkdup.list";
open Metric,">$dOut/metric.list";
my %bam;
my $number=0;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	print Out "$sampleID\t$dOut/$sampleID.mkdup.bam\n";
	print Metric "$sampleID\t$dOut/$sampleID.metric\n";
	print SH "java -Xmx30G -jar /mnt/ilustre/users/dna/.env/bin/picard.jar MarkDuplicates TMP_DIR=$dOut/MKDUP/ MAX_FILE_HANDLES=100 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true I=$bam O=$dOut/$sampleID.mkdup.bam  M=$dOut/$sampleID.metric CREATE_INDEX=TRUE \n";
}
close In;
close SH;
close Out;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=30G --CPU 1 --maxjob $proc  $dShell/05.bam-mkdup.sh";
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
	perl $Script -bam -out -dsh

Usage:
  Options:
  -bam	<file>	input bamlist file
  -out	<out>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	output shell dir

  -h         Help

USAGE
        print $usage;
        exit;
}
