#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fGenome,$fOut,$Gff,$match);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fGenome,
	"o:s"=>\$fOut,
	"g:s"=>\$Gff,
	"f:s"=>\$match
	) or &USAGE;
&USAGE unless ($fGenome and $fOut and $Gff);
open In,$fGenome;
$/=">";
my %seq;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	$seq{$id}=join("\n",@seq);
}
close In;
$/="\n";
open In,$Gff;
open Out,">$fOut";
my %out;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($chrid,$source,$type,$start,$end,undef,undef,undef,$info)=split(/\s+/,$_);
	next if ($type eq "region" || $type eq "exon" || $type eq "CDS");
	my $id="";
	if ($info=~/ID=([^;]*)/) {
		$id=$1;
	}
	if ($info=~/Parent=([^;]*)/) {
		$id=$1;
	}
	next if (exists $out{$id});
	$out{$id}=1;
	print Out ">$id\n";
	print Out substr($seq{$chrid},$start,$end-$start+1),"\n";

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
	reformat genome,rename scaffold name at genome fa file and gff file
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir 

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -o	<str>	output file prefix

  -h         Help

USAGE
        print $usage;
        exit;
}
