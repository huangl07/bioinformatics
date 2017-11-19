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
&USAGE unless ($fGenome and $fOut and $Gff and $match);
my %change;
if (!$match) {
	open In,$match;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$change)=split(/\s+/,$_);
		$change{$id}=$change;
	}
	close In;
}
open In,$fGenome;
open Out,">$fOut.fa"
$/=">";
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($info,@seq)=split(/\n/,$_);
	my $id=(split(/\s+/,$info))[0];
	if (!exists $change{$id}) {
	}else{
		$n++;
		$change{$id}="sca$n";
		print ">$change{$id}";
		print join("\n",@seq),"\n";
	}
	print ">$change{$id}\n";
	print join("\n",@seq),"\n";
}
close In;
close Out;
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
  -f	<file>	chromosome change file

  -h         Help

USAGE
        print $usage;
        exit;
}