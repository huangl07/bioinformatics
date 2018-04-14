#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$mark,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($vcf and $mark );
#$dOut||="./";
my %mark;
open In,$mark;
open Type,">$fOut";
#open Vcf,">$fOut/3-13.xls";
while (<In>){
	next if ($_ eq /#/|| /^$/|| /sca/);
	chomp;
	my ($id,$type,@undi)=split(/\t/,$_,3);
	my $lg=(split(/\D+/,(split(/\_/,$id))[0]))[-1];
	$mark{$lg}{$id}{type}=$type;
}
close In;
open In,$vcf;
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/|| /sca/);
	my ($chr,undef,$id,$ref,$alt,@undi)=split(/\t/,$_,7);#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	my $lg=(split(/\D+/,$chr))[-1];
	$mark{$lg}{$id}{ref}=length($ref);
	$mark{$lg}{$id}{alt}=length($alt);
}
close In;
my ($indel,$snp);
my($Saabb,$Sabcd,$Sabcc,$Sccab,$Sefeg,$Shkhk,$Snnnp,$Slmll,$Iaabb,$Iabcd,$Iabcc,$Iccab,$Iefeg,$Ihkhk,$Innnp,$Ilmll);
#print Vcf "#LG\tSNP Number\tInDel Number\n";
print Type "#TYPE\tSNP Number\tInDel Number\n";
foreach my $lg (sort{$a<=>$b}keys %mark){
	foreach my $id (keys %{$mark{$lg}}){
		if ($mark{$lg}{$id}{ref}==1 and $mark{$lg}{$id}{alt}==1){
			$mark{$lg}{snp}++ ;
			$Saabb++ if($mark{$lg}{$id}{type}=~/aaxbb/);
			$Sabcd++ if($mark{$lg}{$id}{type}=~/abxcd/);
			$Sabcc++ if($mark{$lg}{$id}{type}=~/abxcc/);
			$Sccab++ if($mark{$lg}{$id}{type}=~/ccxab/);
			$Sefeg++ if($mark{$lg}{$id}{type}=~/efxeg/);
			$Shkhk++ if($mark{$lg}{$id}{type}=~/hkxhk/);
			$Snnnp++ if($mark{$lg}{$id}{type}=~/nnxnp/);
			$Slmll++ if($mark{$lg}{$id}{type}=~/lmxll/);
		}else{
			$mark{$lg}{indel}++ ;
			$Iaabb++ if($mark{$lg}{$id}{type}=~/aaxbb/);
			$Iabcd++ if($mark{$lg}{$id}{type}=~/abxcd/);
			$Iabcc++ if($mark{$lg}{$id}{type}=~/abxcc/);
			$Iccab++ if($mark{$lg}{$id}{type}=~/ccxab/);
			$Iefeg++ if($mark{$lg}{$id}{type}=~/efxeg/);
			$Ihkhk++ if($mark{$lg}{$id}{type}=~/efxeg/);
			$Innnp++ if($mark{$lg}{$id}{type}=~/nnxnp/);
			$Ilmll++ if($mark{$lg}{$id}{type}=~/lmxll/);
		}
	}
	#print Vcf "$lg\t$mark{$lg}{snp}\t$mark{$lg}{indel}\n";
}
print Type "aaxbb\t$Saabb\t$Iaabb\n";
print Type "abxcd\t$Sabcd\t$Iabcd\nabxcc\t$Sabcc\t$Iabcc\nccxab\t$Sccab\t$Iccab\nefxeg\t$Sefeg\t$Iefeg\nhkxhk\t$Shkhk\t$Ihkhk\nnnxnp\t$Snnnp\t$Innnp\nlmxll\t$Slmll\t$Ilmll\n";
#close Vcf;
close Type;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input pop.vcf file
  -mark	<file>  input marker file 
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}
