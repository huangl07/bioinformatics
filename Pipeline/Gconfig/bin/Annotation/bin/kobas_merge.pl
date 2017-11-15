#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$dIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($dIn and $fOut);
my @kobas=glob("$dIn/*.kobas");
open Out,">$fOut";
print Out "#query\tKID\tKID_anno\tKoID\tKoanno\n";
foreach my $kobas (@kobas) {
	$/="///";
	open In,$kobas;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my @line=split(/\n/,$_);
		next if (scalar @line ==2);
		my (undef,$queid,$KID,$pathway)=split(/\n/,$_,4);
		$queid=(split(/\t/,$queid))[-1];
		$KID=(split(/\t/,$KID,2))[-1];
		my @Path=split(/\n/,$pathway);
		my @koid;
		my @anno;
		foreach my $path (@Path) {
				my @info=split(/\t/,$path);
				push @koid,$info[3];
				push @anno,$info[1];
		}
		print Out $queid,"\t",$KID,"\t",join(":",@koid),"\t",join(":",@anno),"\n";
	}
	close In;
}
close Out;



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
  -i	<dir>	input *.kobas file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
