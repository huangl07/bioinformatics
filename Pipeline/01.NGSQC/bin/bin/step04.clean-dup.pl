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
open SH,">$dShell/step04.clean-dup.sh";
open Dup,">$dOut/dup.list";
mkdir "$dOut/fig" if (!-d "$dOut/fig");
my %lane;
my $maxmem=0;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2,$e1,$e2)=split(/\s+/,$_);
	print SH "$Bin/bin/ngsduplicate -1 $fq1 -2 $fq2 -q 33 -s $dOut/$sample.dup\n";
	my $mem1=`du $fq1`;
	my $mem2=`du $fq2`;
	$mem1=~s/\D//g;
	$mem2=~s/\D//g;
	if ($mem1+$mem2 >$maxmem) {
		$maxmem=$mem1+$mem2;
	}
	print Dup "$sample\t$dOut/$sample.dup\n";
}
close In;
close SH;
close Dup;
my $mem=join("",int($maxmem/1000000)+1)."G";
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dShell/step04.clean-dup.sh --Resource mem=3G --CPU 1 --Maxjob 19";
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
	-fqlist	<file>	input fqlist file 
	-out	<dir>	output dir
	-dsh	<dir>	output shell	
USAGE
	print $usage;
	exit;
}
