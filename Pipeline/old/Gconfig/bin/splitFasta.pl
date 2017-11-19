#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Split,$Key);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"n:s"=>\$Split,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $Key and $Split);
$Split||=2000;
open In,$fIn;
mkdir $dOut if (!-d $dOut);
$dOut=Cwd::abs_path($dOut);
my %handsh;
$/=">";
my $nseq=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	$nseq++;
	my $filehand=$nseq % $Split;
	if (!exists $handsh{$filehand}) {
		open $handsh{$filehand},">$dOut/$Key.$filehand.fa";
	}
	my ($id,@line)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	print {$handsh{$filehand}} ">$id\n",join("\n",@line),"\n";
}
close In;
foreach my $key (sort keys %handsh) {
	close $handsh{$key};
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	this Script will split fasta file into n files

	eg:
	perl $Script -i demo.fa -d Nr -o ./ -k demo

Usage:
  Options:
  -i	<file>	input fa file name
  -o	<dir>	output dir 
  -k	<str>	output keys of filename
  -n	<num>	split file number
  -h         Help

USAGE
        print $usage;
        exit;
}
