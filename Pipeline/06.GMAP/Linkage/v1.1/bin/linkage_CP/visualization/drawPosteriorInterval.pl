#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$dOut,$log);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($fIn and $fKey);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------
my (@posterior_position) = ();

#-------------------------------------------------------------------
# Get Data
#-------------------------------------------------------------------
open (IN,"$fIn") or die $!;
my $order = 0;
while (<IN>) {
	chomp;
	next if (/^$/ || /^;/) ;
	my ($marker, $type, $position, @plausible) = split /\s+/, $_;
	push @{$posterior_position[$order++]}, @plausible;
}
close (IN) ;
#-------------------------------------------------------------------
# Process
#-------------------------------------------------------------------
open (OUT,">$dOut/$fKey.list") or die $!;

my $xMax = @posterior_position;
my $yMax = $xMax;

#### Distribution 

print OUT "Type:point\n";
print OUT "PointSize:6\n";
print OUT "FontSize:32\n";
print OUT "Width:800\n";
print OUT "Height:800\n";
print OUT "WholeScale:1\n";
print OUT "XStart:0\n";
print OUT "XStep:5\n";
print OUT "XEnd:$xMax\n";
print OUT "YStart:0\n";
print OUT "YEnd:$yMax\n";
print OUT "YStep:5\n";
print OUT "XUnit:0.8\n";
print OUT "MarkPos:rt\n";
print OUT "MarkScale:1\n";
print OUT "Note:Distribution of posterior intervals\n";
print OUT "X:'best' map position\n";
print OUT "Y:plausible map position\n";
print OUT "\n\n";

print OUT "Color:blue\n";
print OUT "NoConnect:1\n";

for (my $i=0;$i<@posterior_position ;$i++) {
	for (my $j=0;$j<@posterior_position ;$j++) {
		if ($posterior_position[$i][$j]) {
			print OUT "$i:$j\n";
		}
	}
}
close (OUT) ;


my $pwd = `pwd`; chomp $pwd;
chdir $dOut;

`perl $Bin/distributing_svg.pl $fKey.list $fKey.list.svg`;
### convert to png 
`$Bin/svg2xxx_release/svg2xxx $fKey.list.svg`;

chdir $pwd;

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
  -i	<file>	Input file, ppm format, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -h		Help

USAGE
	print $usage;
	exit;
}