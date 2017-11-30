#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fPwd,$fKey,$dOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fIn and $fKey );
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$fKey="$dOut/$fKey";

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

open (MST , $fIn) or die $!;
open (MAT , ">$fKey.heatMap") or die $!;

my $flag = 0 ;
while (<MST>) {
	chomp;
	if (/==========/) {
		$flag++ ;
		if ($flag ==2 ) {
			print MAT "R";
		}
		next;
	}

	if ($flag == 2 ) {
		print MAT $_ ,"\n" ; 
	}else{
		next ;
	}
}
close (MST) ;
close (MAT) ;

`$Bin/matrix2png -data $fKey.heatMap -rg -size 8:8 -con 0.75 -d -mincolor yellow -maxcolor blue  -midcolor red  >$fKey.heatMap.png`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:

Usage:
  Options:
  -i	<file>	mst.o file, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}