#!/usr/bin/perl -w
# Writer:       liangsh <liangsh@biomarker.com.cn>
# Last Modified:  2015/10/9.

use strict;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#require bmkPerlBase;
use Cwd qw(abs_path);
use lib "/share/nas33/liangsh/testing/Statistics-RankCorrelation-0.1204/lib/Statistics/";
use RankCorrelation;
my $version = "1.0";
# ------------------------------------------------------------------
# main pipeline
# ------------------------------------------------------------------


my ($mapfile,$locfile,$posifile,$fout);
GetOptions(
			"help|?"=>\&USAGE,
			"m:s"=>\$mapfile,
			"l:s"=>\$locfile,
			"p:s"=>\$posifile,
			"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($mapfile && $locfile && $posifile && $fout);

open (MAP,">$fout.ref.map") or die $!;
open (LOC,">$fout.ref.order.loc") or die $!;

###########get posi Information

open (IN,"<$posifile") or die $!;
my %posi;
while (<IN>) {
	chomp;
	my ($ID,$chr,$str,$end) = split;
	$posi{$ID} = $str;
}
close IN;

###########get initial Map

open (IN,"<$mapfile") or die $!;
my (@marker,@dis,@group);
while (<IN>) {
	chomp;
	next if (/group/);
	my ($mar,$dis) = split;
	push @marker,$mar;
	push @dis,$dis;
	push @group,$mar,$dis;
}
close IN;

#############get Locfile information

open (IN,"<$locfile") or die $!;
my %loc;
while (<IN>) {
	chomp;
	if (! /Marker/) {
		print LOC "$_\n";
		next;
	}
	my ($marker,@ph) = split;
	$loc{$marker} = \@ph;
}
close IN;

#############################output the map and loc with right order
my @ref;
for (my $i=0; $i<@marker; $i++) {
	if (! $posi{$marker[$i]}) {
		die "The posifile is wrong!\n";
	}
	push @ref,$posi{$marker[$i]};
}
my $c = Statistics::RankCorrelation->new(\@dis, \@ref);
my $n = $c->spearman;
#print "$n\n";
if ($n>0) {
	print MAP "group 0\n";
	for (my $i=0; $i<@marker; $i++) {
		print MAP "$marker[$i]\t$dis[$i]\n";
		print LOC $marker[$i],"\t",join("\t",@{$loc{$marker[$i]}}),"\n";
	}
}

if ($n<0) {
	print MAP "group 0\n";
	printf MAP "%s\t%.3f\n",$group[-2],0.000;
	print LOC $group[-2],"\t" , join("\t",@{$loc{$group[-2]}}) , "\n" ;
	my $sum;
	for (my $i=2; $i<@group; $i+=2) {
			$sum += ($group[-$i+1]-$group[-$i-1]);
			$sum =~ s/^\s+//;
			printf MAP "%s\t%.3f\n",$group[-$i-2],$sum;
			print LOC $group[-$i-2],"\t" , join("\t",@{$loc{$group[-$i-2]}}),"\n" ;
		}
}

close MAP;
close LOC;

unlink ("$mapfile") or die $!;
unlink ("$locfile") or die $!;
rename ("$fout.ref.map","$mapfile") or die $!;
rename ("$fout.ref.order.loc","$locfile") or die $!;







# ------------------------------------------------------------------
# print usage information and exit
# ------------------------------------------------------------------
sub USAGE {
	my $usage=<<"USAGE";
Version:	v$version
Writer:		liangsh <liangsh\@biomarker.com.cn>
Program Date:	2015/01/01.
Description:	this program is code for order the mapfile and locfile by the reference
Usage:
  Options:
  Forced parameters:
  -m	<str>	the name of mapfile
  -l	<str>	the name of locfile
  -p	<str>	the name of posifile
  -o	<str>	the outputdir
  Optional parameters:
  -h		<none>	Help
USAGE
	print $usage;
	exit(1);
}




