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
open SH,">$dShell/04.bam-mkdup.sh";
open In,$bamlist;
open Out,">$dOut/bam.list";
open Metric,">$dOut/metric.list";
my %bam;
my $number=0;
my $sentieon="/mnt/ilustre/users/dna/.env/sentieon/sentieon/sentieon-genomics-201711.05/bin/sentieon";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam,$type)=split(/\s+/,$_);
	#print Out "$sampleID\t$dOut/$sampleID.mkdup.bam\t$dOut/$sampleID.recal_data.txt\n";
	print Out "$sampleID\t$dOut/$sampleID.mkdup.bam\n";
	print Metric "$sampleID\t$dOut/$sampleID.metric.txt\n";
	if ($type ne "WGS") {
		print SH "ln -s $bam $dOut/$sampleID.mkdup.bam \n";
	}else{
		if (!-f $bam) {
			die "check $bam!";
		}
		print SH "$sentieon driver -t 1 -i $bam --algo LocusCollector --fun score_info $dOut/$sampleID.score.txt && ";
		print SH "$sentieon driver -t 8 -i $bam --algo Dedup --rmdup --score_info $dOut/$sampleID.score.txt --metrics $dOut/$sampleID.metric.txt $dOut/$sampleID.mkdup.bam \n";
	}
	#print SH "/mnt/ilustre/users/dna/.env/sentieon/sentieon/sentieon-genomics-201711.05/bin/sentieon driver -t 1 -r $ref -i $dOut/$sampleID.mkdup.bam --algo QualCal  $dOut/$sampleID.recal_data.txt\n";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=30G --CPU 1 --maxjob $proc $dShell/04.bam-mkdup.sh";
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
