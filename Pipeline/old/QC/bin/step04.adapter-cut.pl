#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fqlist,$dOut,$dShell);
GetOptions(
				"help|?" =>\&USAGE,
				"fqlist:s"=>\$fqlist,,
				"out:s"=>\$dOut,
				"dsh:s"=>\$dShell,
				) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);
open In,$fqlist;
open SH,">$dShell/step04.adapter-cut.sh";
mkdir "$dOut/fig" if (!-d "$dOut/fig");
open Out,">$dOut/fastq.list";
open Stat,">$dOut/ada.list";
my %lane;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2,$e1,$e2)=split(/\s+/,$_);
	print SH "/mnt/ilustre/users/dna/.env//bin/SeqPrep  -f $fq1 -r $fq2 -1 $dOut/$sample.rm_1.fastq.gz -2 $dOut/$sample.rm_2.fastq.gz ";
	print SH "-q 20 -L 25 -B AGATCGGAAGAGCGTCGTGT -A AGATCGGAAGAGCACACGTC 2>$dOut/$sample.rm_adapter.log \n";
	print Out "$sample $dOut/$sample.rm_1.fastq.gz $dOut/$sample.rm_2.fastq.gz\n";
	print Stat "$sample\t$dOut/$sample.rm_adapter.log\n";
}
close In;
close SH;
close Out;
close Stat;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dShell/step04.adapter-cut.sh --Resource mem=20G";
`$job`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-fqlist	<file>	input file 
	-out	<dir>	output dir
	-dsh	<dir>	output work sh
USAGE
	print $usage;
	exit;
}
