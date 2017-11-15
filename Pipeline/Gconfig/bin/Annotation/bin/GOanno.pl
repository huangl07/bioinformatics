#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut";
print Out "#geneid\tGO_id\tGO_trans\n";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my @info=split(/\t/,$_);
	my @GOid;
	my @detail;
	while ($info[-1] =~ m/\[(GO:\d+) "([^\"\']*)" evidence=\w+\]/g) {
		push @GOid,$1;
		push @detail,$2;
	}
	if (scalar @GOid == 0 || scalar @detail == 0) {
		print Out $info[0],"--\t--\n";
		next;
	}
	print Out $info[0],"\t",join(":",@GOid),"\t",join(":",@detail),"\n";
}
close Out;
close In;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
