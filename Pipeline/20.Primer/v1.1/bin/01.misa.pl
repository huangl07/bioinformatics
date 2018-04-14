#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
		"help|?" =>\&USAGE,
		"i|input:s"=>\$vcf,
		"o|output:s"=>\$out,
		"t|type:s"=>\$type,
                        ) or &USAGE;
&USAGE unless ($vcf and $out and $type);
$type||="SNP";
my ($totalnum,$sum_chr);
my $Chr="";
open In,$vcf;
open Out,">$out.$type.misa";
print Out "#chr:ID\t$type.nr\ttotalnumber\ttype\tsnp\talt\tsize\tstart\tend\tnstart\tnend\n";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/);
	my ($chr,$pos,$ids,$ref,$alt,undef)=split(/\s+/,$_);
	$ids=~s/\_/\:/;
	my $chrtype=(split(/\d/,$chr))[0];
	my $chrnr=(split(/\D/,$chr))[-1];
	my $size=length($ref);
	my $start=$pos;
	my $end=$pos+$size - 1;
	$totalnum++;					## if SNP/INDEL,need seprately sum;
	my($nstart,$nend);
	if ($type=~/SNP/i){     
############################################################################### change length
		$nstart=$pos - 500 ;
		$nend=$pos+  500 ;
	}else{
		$nstart=$pos - 150 ;
		$nend=$pos+$size+149;
	}
	if ($Chr ne $chr){
		$sum_chr=1;
		$Chr=$chr;
	}elsif($Chr eq $chr){
		$sum_chr++;
	}
	print Out join("\t",$ids,$sum_chr,$totalnum,$type,$ref,$alt,$size,$start,$end,$nstart,$nend),"\n";
}	
close In;
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
		-input	<file>  input vatiant vcf file 
		-output	<dir>  output misa dir
		-type	<str>	marker type
USAGE
        print $usage;
        exit;
}
