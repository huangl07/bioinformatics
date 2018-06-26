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
open In1,$fIn1;
open In2,$fIn2;
open Out,">",$fOut;
my (%hash1,%hash2,%hash3);
while (<In1>){
	chomp;
	next if (/^$/ or "" or /#/);
	my @line=split(/\t/);
	my $name=join("",$line[0],$line[1]);
	$hash1{$name}=1;
	$hash2{$name}=$_;
}
close In1;
while(<In2>){
	chomp;
	next if (/^$/ or "" or /#/);
	my @line=split(/\t/);
	if ($line[0] =~ /ref/){
		my $name=join("",$line[0],$line[1]);
		$hash2{$name}=$_;
		$hash3{$name}=$_;
	}
	else{
		my $name=join("","seq",$line[0]);
		$hash2{$name}=$_;
		$hash3{$name}=$_;
	}
}
close In2;
foreach my $key (sort {$a cmp $b} keys %hash2){
	print Out"$hash2{$key}";
	if (exists $hash1{$key} and exists $hash3{$key}){
		print Out"\t1\/1\t1\/1\n";
	}
	if (!exists $hash1{$key} and exists $hash3{$key}){
		print Out"\t0\/0\t1\/1\n";
	}
	if (exists $hash1{$key} and !exists $hash3{$key}){
		print Out"\t1\/1\t0\/0\n";
	}
	if (!exists $hash1{$key} and !exists $hash3{$key}){
		print Out"\t0\/0\t0\/0\n";
	}
}
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
  -i	<file>	misa.out1 file name
  -k	<file>	misa.out1 file name
  -o	<file>	out file name
  -h         Help

USAGE
        print $usage;
        exit;
}
