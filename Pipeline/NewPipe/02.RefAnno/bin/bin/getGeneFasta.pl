#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$ref,
	"o:s"=>\$fOut,
	) or &USAGE;
&USAGE unless ($ref and $out and $Gff);
my %enzyme=(
	"EcoRI"=>"GAATTC",
	"MseI"=>"TTAA",
	"PstI"=>"CTGCAG",
	"TaqaI"=>"TCGA",
);
my @enzyme=keys %enzyme;
open In,$ref;
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	my $seq=join("",@seq);
	
}
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
