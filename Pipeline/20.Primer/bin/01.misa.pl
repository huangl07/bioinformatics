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
		"input:s"=>\$vcf,
		"type:s"=>\$type,
		"output:s"=>\$out,
                        ) or &USAGE;
&USAGE unless ($vcf  and $out);
#mkdir $out if (!-d $out);
my (@unif,$start,$end,$ids,$chr,$pos,$size,$nstart,$nend);
my $num=1;
my $number=0;
my $chrnum=1;
my $scanum=1;
my $newnum=1;
open In,$vcf;
open Out,">$out.$type.misa";
print Out "#chr:ID\t$type.nr\ttotalnumber\ttype\tsnp\talt\tsize\tstart\tend\tnstart\tnend\n";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/);
	my ($chr,$pos,$ids,$ref,$alt,@unif)=split(/\s+/,$_);
	$ids=~s/\_/\:/;
	my $chrtype=(split(/\d/,$chr))[0];
	my $chrnr=(split(/\D/,$chr))[-1];
	$size=length($ref);
	$start=$pos;
	$end=$pos+$size - 1;
	$number++ ;
	if ($type=~/snp/){
		$nstart=$pos - 500 ;
		$nend=$pos+  500 ;
	}else{
		$nstart=$pos - 150 ;
		$nend=$pos+$size+149;
	}
	if ($chrtype=~/chr/){
		if ($chrnr eq $chrnum){
			print Out join("\t",$ids,$num,$number,$type,$ref,$alt,$size,$start,$end,$nstart,$nend),"\n";
			$num++;
		}else{
			$num=1;
			print Out join("\t",$ids,$num,$number,$type,$ref,$alt,$size,$start,$end,$nstart,$nend),"\n";
			$num++;
			$chrnum= $chrnum +1;
		}
	}
	if ($chrtype=~/sca/){
		if ($chrnr eq $scanum){
			print Out join("\t",$ids,$newnum,$number,$type,$ref,$alt,$size,$start,$end,$nstart,$nend),"\n";
			$newnum++;
		}else{
			$newnum=1;
			print Out join("\t",$ids,$newnum,$number,$type,$ref,$alt,$size,$start,$end,$nstart,$nend),"\n";
			$newnum++;
			$scanum=$scanum +1;
		}
	}
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
