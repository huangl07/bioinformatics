#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($type,$vcf,$mark,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"type:s"=>\$type,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($vcf and $mark );
#$fOut||="./";
my %mark;
$type||="chr";
open In,$mark;
open Vcf,">$fOut";
while (<In>){
	next if ($_=~/^#/|| /^$/|| /group/);
	chomp;
	my ($id,@undi)=split(/\t/,$_);
	my $lg=(split(/\D+/,(split(/\_/,$id))[0]))[-1];
	my $pos=(split(/\_/,$id))[1];
	$mark{$lg}{$id}{pos}=$pos;
}
close In;
open In,$vcf;
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/|| /sca/);
	my ($chr,$pos,$id,$ref,$alt,@undi)=split(/\t/,$_,7);#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	#my $lg=(split(/\D+/,$chr))[-1];
	if($chr=~/$type/){
		my $lg=(split(/\D+/,$chr))[-1];
		$mark{$lg}{$id}{ref}=length($ref);
		$mark{$lg}{$id}{alt}=length($alt);
		$mark{$lg}{$id}{vcpos}=$pos;
	}else{
		next;
	}
}
close In;
my ($indel,$snp);
print Vcf "#LG\tSNP Number\tInDel Number\n";
foreach my $lg (sort{$a<=>$b}keys %mark){
	foreach my $id (keys %{$mark{$lg}}){
		if($mark{$lg}{$id}{pos} eq $mark{$lg}{$id}{vcpos}){
			if ($mark{$lg}{$id}{ref}==1 and $mark{$lg}{$id}{alt}==1){
				$mark{$lg}{snp}++ ;
			}else{
				$mark{$lg}{indel}++ ;
			}
		}
	}
	print Vcf "$lg\t$mark{$lg}{snp}\t$mark{$lg}{indel}\n";
}
close Vcf;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
  -vcf	<file>	input pop.vcf file
  -type	<stri>	chr	or sca	
  -mark	<file>  input marker file 
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}
