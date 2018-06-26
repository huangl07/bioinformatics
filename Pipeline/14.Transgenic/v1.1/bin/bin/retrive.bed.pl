#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($sam,$id,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "sam:s"=>\$sam,
    "id:s"=>\$id,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($sam and $id and $out);
mkdir $out if (!-d $out);
open IN,"$sam";
#open IN,"$bam";
open OUT,">$out/$id.bed";

chomp(my @data=<IN>);
for(my $i=0;$i<=$#data;$i++){
    my ($fqid1,$flag1,$sca1_1,$pos1_1,$a1,$b1,$sca1_2,$pos1_2,undef,$seq1,$qual1,@others1)=split(/\t/,$data[$i]);
    my ($fqid2,$flag2,$sca2_1,$pos2_1,$a2,$b2,$sca2_2,$pos2_2,undef,$seq2,$qual2,@others2)=split(/\t/,$data[$i+1]);
    if($sca1_1 eq "Chrinsert" and $sca2_1=~/sca/){
        my $start=$pos2_1-500;
        my $end=$pos2_1+500;
        #open OUT,">$out/$id.$sca1_1.$pos1_1.$sca2_1.$pos2_1.bed";
        print OUT "$fqid1\t$sca1_1\t$sca2_1\t$start\t$end\n";
        #close OUT;
        }elsif($sca1_1=~/sca/ and $sca2_1 eq "Chrinsert"){
        my $start=$pos1_1-500;
        my $end=$pos1_1+500;
        #open OUT,">$out/$id.$sca1_1.$pos1_1.$sca2_1.$pos2_1.bed";
        print OUT "$fqid1\t$sca1_1\t$sca2_1\t$start\t$end\n";
        #close OUT;
        }else{
        #open OUT,">$out/$id.Chrinsert.bed";
        print OUT "$fqid1\t$sca1_1\t$sca2_1\t$pos1_1\t$pos2_1\n";
        #close OUT;
        }
    }
close OUT;




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
  -sam    <file>  input bam file
  -id   bam id 
  -out	<dir>	out dir     
  -h         Help

USAGE
        print $usage;
        exit;
}
