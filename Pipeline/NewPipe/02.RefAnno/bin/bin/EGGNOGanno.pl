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
		my %EggNOG=(
			"J"=>"Translation, ribosomal structure and biogenesis",
			"A"=>"RNA processing and modification",
			"K"=>"Transcription",
			"L"=>"Replication, recombination and repair",
			"B"=>"Chromatin structure and dynamics",
			"D"=>"Cell cycle control, cell division, chromosome partitioning",
			"Y"=>"Nuclear structure",
			"V"=>"Defense mechanisms",
			"T"=>"Signal transduction mechanisms",
			"M"=>"Cell wall/membrane/envelope biogenesis",
			"N"=>"Cell motility",
			"Z"=>"Cytoskeleton",
			"W"=>"Extracellular structures",
			"U"=>"Intracellular trafficking, secretion, and vesicular transport",
			"O"=>"Posttranslational modification, protein turnover, chaperones",
			"C"=>"Energy production and conversion",
			"G"=>"Carbohydrate transport and metabolism",
			"E"=>"Amino acid transport and metabolism",
			"F"=>"Nucleotide transport and metabolism",
			"H"=>"Coenzyme transport and metabolism",
			"I"=>"Lipid transport and metabolism",
			"P"=>"Inorganic ion transport and metabolism",
			"Q"=>"Secondary metabolites biosynthesis, transport and catabolism",
			"R"=>"General function prediction only",
			"S"=>"Function unknown",
		);

open In,$fIn;
open Out,">$fOut";
print Out "#geneid\tGO_id\tGO_trans\n";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my @info=split(/\t/,$_);
	print Out $info[0],"\t",$info[-2],"\t",$EggNOG{$info[-2]},"\n";
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
