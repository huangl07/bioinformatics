#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bam,$id,$bed,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "bam:s"=>\$bam,
    "id:s"=>\$id,
    "bed:s"=>\$bed,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bam and $bed and $id and $out);
mkdir $out if (!-d $out);
my %hash;
my %stat;
open BED,"$bed";
while(<BED>){
    $_=~s/[\n\r]//g;
    my ($ids,$chr1,$chr2,$start,$end)=split/\t/,$_;
    if($chr1 eq "Chrinsert"){
        $hash{$chr2}{"$start\t$end"}=1
        }
    if($chr2 eq "Chrinsert"){
        $hash{$chr1}{"$start\t$end"}=1
            }
    }
close BED;

open IN,"samtools view $bam|";

while(<IN>){
    $_=~s/[\n\r]//g;
    my ($fqid,$flag,$sca1,$pos1,$a,$b,$sca2,$pos2,undef,$seq,$qual,@others)=split(/\t/,$_);
    foreach my $pos (sort keys %{$hash{$sca1}}){
        my ($s,$e)=split/\t/,$pos;
        open OUT,">>$out/$id.$sca1.$s.$e.id";
        if($pos1 > $s and $pos1 < $e){
            $stat{$sca1}{$pos}{seq}.=$seq;
            $stat{$sca1}{$pos}{qual}.=$qual;
            $stat{$sca1}{$pos}{num}++;
            print OUT "$fqid\t$sca1\t$pos1\t$sca2\t$pos2\n";
            }
        }
    }
close IN;
close OUT;
open ST,">$out/$id.stat.xls";
foreach my $chr (sort keys %stat){
    foreach my $pos (sort keys %{$stat{$chr}}){
        my ($start,$end)=split/\t/,$pos;
        my $length=length($stat{$chr}{$pos}{seq});
        my ($q20,$q30)=qual_stat($stat{$chr}{$pos}{qual});
        my $q30s=sprintf("%4.f",$q30/$length);
        print ST "$chr\t$start\t$end\t$stat{$chr}{$pos}{num}\t$length\t$q30s\n";
        }
    }
sub qual_stat{
    my @qstr=split//,$_;
    my $q20;
    my $q30;
    for(my $i=0;$i<=$#qstr;$i++){
        my $qual = ord($qstr[$i])-33;
        if($qual >=30){
            $q30++;
            $q20++;
       }
        if($qual >=20){
            $q20++;
       }
     }
    return $q20,$q30
}
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
  -bed  <dir> dir;bed file from retrive bed
  -out	<dir>	out dir     
  -h         Help

USAGE
        print $usage;
        exit;
}
