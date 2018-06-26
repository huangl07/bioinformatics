#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bam,$id,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "bam:s"=>\$bam,
    "id:s"=>\$id,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bam and $id and $out);
mkdir $out if (!-d $out);
open IN,"samtools view $bam|";
open OUT,">$out/$id.insert.sam";
open BED,">$out/$id.bed";
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($fqid,$flag,$sca1,$pos1,$a,$b,$sca2,$pos2,$c,$seq,$qual,@others)=split/\t/,$_;
    if($sca1 eq "Chrinsert" or $sca2 eq "Chrinsert"){
        my $start=$pos1-500;
        my $end=$pos1+500;
        print BED "$fqid\t$sca1\t$sca2\t$start\t$end\n";
        print OUT join("\t",$fqid,$flag,$sca1,$pos1,$a,$b,$sca2,$pos2,$c,$seq,$qual,@others),"\n";
        }
}
close IN;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -g -i -o 

Usage:
  Options:
  -bam    <file>  input bam file
  -id   bam id 
  -out	<dir>	out dir     
  -h         Help

USAGE
        print $usage;
        exit;
}
