#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fIn2,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn1,
	"k:s"=>\$fIn2,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn1 );
open Read,$fIn1;
my %noMatchChr;
while(<Read>){
	chomp;
	next if (/^$/ or "");
	my ($chr_pos,@long)=split(/\s+/);
	$chr_pos=~s/\#//;
	my ($chr,$pos)=split(/\:/,$chr_pos);
	$noMatchChr{$chr}=1;
}
close Read;
$/=">";
open Write,$fIn2;
open Out,">",$fOut;
while(<Write>){
	chomp;
	next if (/^$/ or "");
	my ($chr,@line)=split(/\n/);
	print "$chr\n";die;
	if(exists $noMatchChr{$chr}){
		print Out"\>$_";
	}
}
close Write;
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
