#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$dOut,$dShell);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut = ABSOLUTE_DIR($dOut);
mkdir $dShell if(!-d $dShell);
$dShell = ABSOLUTE_DIR($dShell);
open SH,">$dShell/step01.ustacks.sh";
open Out,">$dOut/ustacks.list";
open In,$fqlist;
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$fq1)=split(/\s+/,$_);
	print SH "ustacks -f $fq1 -o $dOut/ -i $n --deleverage -M 6 -m 2 -p 8 --gapped\n";
	print Out "$sample\t$dOut/$sample\n";
	$n++;
}
close In;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=10G --CPU 8 --Nodes 1 $dShell/step01.ustacks.sh";
print "$job\n";
`$job`;
print "$job\tdone!\n";
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
  -fqlist	<file>	input fqlist name
  -out	<dir>	split windows sh
  -dsh	<dir>	output work sh	
  -h         Help

USAGE
        print $usage;
        exit;
}
