#!/usr/bin/perl -w
use FindBin qw($Bin $Script);
BEGIN{ #add absolute path of required moduals to array \@INC
	push @INC,"$Bin/blib/lib";
}
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Statistics::Distributions;
use Statistics::Regression;
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
				"p:s"=>\$fPwd,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $fPwd);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$fKey="$dOut/$fKey";
#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (@marker,%pwd) = ();

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------
open (M,"$fIn") or die $!;

while (<M>) {
	s/\r//g;
	chomp;
	next if (/^$/ || /^group/ || /^;/) ;
	next unless (/^(\S+)\s+\S+/) ;

	push @marker,$1;

}
close (M) ;

#print Dumper @marker;die;


open (P,"$fPwd") or die $!;

while (<P>) {
	chomp;
	next if (/^$/ || /^group/ || /^;/) ;
	next unless (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) ;

	$pwd{$1}{$2}{'r'}=$3;
	$pwd{$2}{$1}{'r'}=$3;
	$pwd{$1}{$2}{'lod'}=$4;
	$pwd{$2}{$1}{'lod'}=$4;

}

close (P) ;

#print Dumper %pwd;die;

#-------------------------------------------------------------------
# Process
#-------------------------------------------------------------------

open (R,">$fKey.r.heatMap") or die $!;
open (LOD,">$fKey.lod.heatMap") or die $!;

print R join("\t",("R",@marker)),"\n";
print LOD join("\t",("LOD",@marker)),"\n";

for (my $i=0;$i<@marker ;$i++) {

	print R $marker[$i];
	print LOD $marker[$i];

	for (my $j=0;$j<@marker ;$j++) {

		if ($marker[$i] eq $marker[$j]) {

			$pwd{$marker[$i]}{$marker[$j]}{'r'}||=0;
			$pwd{$marker[$i]}{$marker[$j]}{'lod'}||=30;

		}
#		else{
#
#			$pwd{$marker[$i]}{$marker[$j]}{'r'}||=0.5;
#			$pwd{$marker[$i]}{$marker[$j]}{'lod'}||=0;
#		}

		if (exists $pwd{$marker[$i]}{$marker[$j]}) {

			print R "\t",sprintf("%.3f",$pwd{$marker[$i]}{$marker[$j]}{'r'});
			print LOD "\t",sprintf("%.3f",$pwd{$marker[$i]}{$marker[$j]}{'lod'});
		}else{
			
			print R "\t",'-';
			print LOD "\t",'-';
		}

		
	}

	print R "\n";

	print LOD "\n";

}

close (R) ;
close (LOD) ;


`$Bin/matrix2png -data $fKey.r.heatMap -rg -size 8:8 -con 0.75 -d -mincolor yellow -maxcolor blue  -midcolor red  > $fKey.r.heatMap.png`;

`$Bin/matrix2png -data $fKey.lod.heatMap -rg -size 8:8 -con 0.75 -d  -mincolor red -maxcolor blue  -midcolor yellow  > $fKey.lod.heatMap.png`;



#-------------------------------------------------------------------
# Print
#-------------------------------------------------------------------

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
  -i	<file>	Map file, forced
  -p	<file>	PWD file, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}