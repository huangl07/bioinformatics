#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$gff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$ref,
	"o:s"=>\$out,
	"g:s"=>\$gff,
	) or &USAGE;
&USAGE unless ($ref and $out and $gff);
open In,$ref;
$/=">";
my %seq;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	my $seq=join("",@seq);
	$seq{$id}=$seq;
}
close In;
open In,$gff;
open Out,">$out";
$/="\n";
my $gname;
my $id;
my %out;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	my ($chr,$source,$type,$start,$end,undef,undef,undef,$info)=split(/\t/,$_);
	next if ($type eq "CDS" ||$type eq "exon" || $type eq "region");
	if ($info =~/ID=([^;]*)/) {
		$id=$1;
	}
	my $pos=join("\t",$chr,$start,$end);
	next if (exists $out{$pos});
	print Out ">$id:$chr:$start:$end\n";
	print Out substr($seq{$chr},$start,$end-$start+1),"\n";
	$out{$pos}=1;
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

  -h         Help

USAGE
        print $usage;
        exit;
}
