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
my ($fLG,$fGenotype,$dOut);
GetOptions(
				"help|?" =>\&USAGE,
				"l:s"=>\$fLG,
				"g:s"=>\$fGenotype,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fLG and $fGenotype);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./extractGenotype";
mkdir $dOut unless (-d $dOut) ;
#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%gData) = ();

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------

open (G,"$fGenotype") or die $!;
$/="\n";
while (<G>) {
	chomp;
	next if (/^$/ || /^;/ || /^MarkerID/ || /^name/ || /^popt/ || /^nloc/ || /^nind/);
	my ($marker) = split ;
	$gData{$marker}=$_;

}
close (G) ;

open (LG,"$fLG") or die $!;
$/=">";
while (<LG>) {
	chomp;
	next if (/^$/ || /^;/) ;

	my ($lgHead,@markerLine) = split /\n/,$_;

	my ($key) = $lgHead =~/^(\w+)/;
	my @marker = map {chomp $_;split /\s+/,$_;} @markerLine;

	open (OUT,">$dOut/$key.genotype") or die $!;

	foreach my $marker (@marker) {

		if (exists $gData{$marker}) {

			print OUT $gData{$marker},"\n";

		}else{

			print "error!Info of $marker in $key is not found in file $fGenotype,please check your input!";

			exit(-1);
		}
	}

	close (OUT) ;
}

close (LG) ;

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
			This programme is designed to  extact genotype info according to lg file.

Usage:
  Options:
  -l	<file>	Linkage group file, forced
  -g	<file>	Genotype file ,forced
  -d	<str>	Directory where output file produced,optional,default [./extractGenotype]
  -h		Help

USAGE
	print $usage;
	exit;
}
