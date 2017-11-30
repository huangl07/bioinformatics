#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Fatal qw(open close);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fDetail,$fLoc,$fKey,$dOut,$log);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fDetail,
				"l:s"=>\$fLoc,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fDetail and $fLoc and $fKey);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%loc,%pair_lp) = ();

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------
open (IN,"<",$fLoc);

while (<IN>) {

	chomp;
	next if (/^$/ or /^name/ or /^nloc/ or /^nind/ or /^popt/ or /^;/);

	my ($marker,$type,$lp,@geno) = split;

	if ($lp !~/{..}/) {

		print "loc file should contain the linkage phase info of loci\n";

		exit;
	}

	$lp=~s/{|}//g;
	$type=~s/<|>//g;
	$loc{$marker}{'lp'} = $lp;
	$loc{$marker}{'type'} = $type;
	$loc{$marker}{'geno'} = \@geno;
}

close (IN); 

#-------------------------------------------------------------------
# extract pwd in detail file
#-------------------------------------------------------------------

my @marker = keys %loc;

for (my $i=0;$i<@marker-1 ;$i++) {

	for (my $j=$i+1;$j<@marker ;$j++) {

		my ($plp_1,$mlp_1)= $loc{$marker[$i]}{'lp'} =~/(.)(.)/; 
		my ($plp_2,$mlp_2)= $loc{$marker[$j]}{'lp'} =~/(.)(.)/;
		
		my $plp = ($plp_1 eq $plp_2)?'C':'R';
		my $mlp = ($mlp_1 eq $mlp_2)?'C':'R';

		$pair_lp{$marker[$i]}{$marker[$j]} = $plp.$mlp;
		$pair_lp{$marker[$j]}{$marker[$i]} = $plp.$mlp;
	}
}

my $pair_num = 0;

open (IN,"<",$fDetail);
open (OUT,">","$dOut/$fKey.pwd");

while (<IN>) {

	chomp;
	next if (/^$/ or /^;/);

	my ($marker1,$marker2,$type1,$type2,$integrity,$recCC,$lodCC,$fiCC,$recCR,$lodCR,$fiCR,$recRC,$lodRC,$fiRC,$recRR,$lodRR,$fiRR,$mLOD) = split;

	my ($p11,$p12,$m11,$m12) = $type1=~/<?(.)(.)x(.)(.)>?/;
	my ($p21,$p22,$m21,$m22) = $type2=~/<?(.)(.)x(.)(.)>?/;

	next if (($p11 eq $p12 and $m11 eq $m12) or ($p21 eq $p22 and $m21 eq $m22) or ($p11 eq $p12 and $m21 eq $m22) or ($p21 eq $p22 and $m11 eq $m12)) ;

	if (exists $pair_lp{$marker1}{$marker2}) {

		$pair_num++;

		if ($pair_lp{$marker1}{$marker2} eq 'CC') {

			$recCC = sprintf("%.16e",0.499) if ($recCC >= 0.5) ;
			$mLOD = sprintf("%.16e",0.001) if ($recCC >= 0.49) ;

			print OUT join("\t",($marker1,$marker2,$recCC,$mLOD,$lodCC,$fiCC,$integrity,'CC')),"\n";
		}elsif($pair_lp{$marker1}{$marker2} eq 'CR'){

			
			$recCR = sprintf("%.16e",0.499) if ($recCR >= 0.5) ;
			$mLOD = sprintf("%.16e",0.001) if ($recCR >= 0.49) ;

			print OUT join("\t",($marker1,$marker2,$recCR,$mLOD,$lodCR,$fiCR,$integrity,'CR')),"\n";
		}elsif($pair_lp{$marker1}{$marker2} eq 'RC'){

			$recRC = sprintf("%.16e",0.499) if ($recRC >= 0.5) ;
			$mLOD = sprintf("%.16e",0.001) if ($recRC >= 0.49) ;
			
			print OUT join("\t",($marker1,$marker2,$recRC,$mLOD,$lodRC,$fiRC,$integrity,'RC')),"\n";
		}else{

			$recRR = sprintf("%.16e",0.499) if ($recRR >= 0.5) ;
			$mLOD = sprintf("%.16e",0.001) if ($recRR >= 0.49) ;

			print OUT join("\t",($marker1,$marker2,$recRR,$mLOD,$lodRR,$fiRR,$integrity,'RR')),"\n";
		}
	}
}

close (IN);
close (OUT);

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
  -i	<file>	PWD detail file , forced
  -l	<file>	Loc file with linkage phase, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}