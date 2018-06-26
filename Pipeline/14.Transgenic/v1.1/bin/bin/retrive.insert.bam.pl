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
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bam and $out);
mkdir $out if (!-d $out);
open IN,"samtools view $bam|";
my %filehand;
my %stat;
my %region;
while (<IN>) {
	next if ($_ eq ""||/^$/);
    my ($fqid,$flag,$sca1,$pos1,$a,$b,$sca2,$pos2,$c,$seq,$qual,@others)=split/\t/,$_;
	if ($sca1 eq "Chrinsert" && $sca2 eq "=") {#reads本身存在于插入序列上
		if (!exists $filehand{$Chrinsert}) {
			open $filehand{$Chrinsert},">$out/Chrinsert.fq";
		}
		print {$filehand{$Chrinsert}},"\@$sca1,$pos1\n$seq\n+\n$qual\n";
		$stat{$Chrinsert}{readnum}++;
		$stat{$Chrinsert}{basenum}++;
		$stat{$Chrinsert}{q30}++;
	}else{
		if ($sca1 eq "Chrinsert" || $sca2 eq "Chrinsert") {
			my $pos=$pos1;
			my $chr=$sca1;
			if ($sca1 eq "Chrinsert") {
				$pos=$pos2;
				$chr=$sca2;
			}
			my $region=join("\t",$chr,$pos-500,$pos+500);
			if (scalar keys %region > 0) {
				$region{$region}=1;
			}else{
				my $check=0;
				foreach $region2 (sort keys %region) {
					my ($chr1,$start1,$end1)=split(/\s+/,$region);
					my ($chr2,$start2,$end2)=split(/\s+/,$region2);
					next if ($chr1 ne $chr2);
					if ($start2 > $end1 || $end2 < $start1) {
						$check=0;
						next;
					}else{
						$check=1;
						my $news=$start1;
						if ($start1 > $start2 ) {
							$news=$start2;
						}
						my $newe=$end1;
						if ($end1 < $end2) {
							$newe=$end2;
						}
						delete $region{$region2};
						$region{join("\t",$chr1,$news,$newe)}=1;
					}
				}
				if ($check ==0) {
					$region{$region}=1;
				}
			}
		}
	}
}
close OUT;
print Out "$out/region.bed";
foreach my $region(sort keys %region) {
	print Out $region,"\n";
}
close Out;



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
  -bam	<file>  input bam file
  -out	<dir>	out dir     
  -h         Help

USAGE
        print $usage;
        exit;
}
