#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dIn,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$dIn,
	"o:s"=>\$out,
			) or &USAGE;
&USAGE unless ($dIn and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input trit's qtl.csv file's dir 
  -o	<file>	output qtl's result xls
  -h         Help

USAGE
        print $usage;
        exit;
}

my @file=glob("$dIn/*.qtl.csv");
open OUT,">$out";
#print OUT "Trait\tChr\tPosition(cM)\tLOD\tR2(%)\tStart(cM)\tEnd(cM)\tMarker Number\n";
print OUT "Trait\tChr\tPosition(cM)\tLOD\tR2(%)\tStart(cM)\tEnd(cM)\tMarker ID\tMarker start\tMarker end\n";
my %stat;

for(@file){
    my $trit=(split(/.qtl.csv/,(split(/\//,$_))[-1]))[0];
    open In,"$_";
    while(<In>){
        s/\"//g;
		chomp;
        next if ($_=~/marker/|| /^$/);
        my ($marker,$chr,$pos,$lod,$var,$pm1,$pm2,$start,$end,$mark1,$mark2)=split/\t/,$_;
		#$marker=(split(/\_/,$marker))[-1];
		#$stat{$trit}{$chr}{$marker}{pos}=$pos;
		#$stat{$trit}{$chr}{$marker}{lod}=$lod;
		#$stat{$trit}{$chr}{$marker}{var}=$var;
		#$stat{$trit}{$chr}{$marker}{start}=$start;
		#$stat{$trit}{$chr}{$marker}{end}=$end;
		$stat{$trit}{$chr}{$pos}{marker}=$marker;
		$stat{$trit}{$chr}{$pos}{lod}=$lod;
		$stat{$trit}{$chr}{$pos}{var}=$var;
		$stat{$trit}{$chr}{$pos}{start}=$start;
		$stat{$trit}{$chr}{$pos}{end}=$end;
		$stat{$trit}{$chr}{$pos}{mark1}=$mark1;
		$stat{$trit}{$chr}{$pos}{mark2}=$mark2;
	}
	close In;
}
foreach my $trit (sort keys %stat){
	foreach my $chr (sort{$a<=>$b}keys%{$stat{$trit}}){
		foreach my $pos (sort{$a<=>$b}keys%{$stat{$trit}{$chr}}){
			print OUT join("\t",$trit,$chr,$pos,$stat{$trit}{$chr}{$pos}{lod},$stat{$trit}{$chr}{$pos}{var},$stat{$trit}{$chr}{$pos}{start},$stat{$trit}{$chr}{$pos}{end},$stat{$trit}{$chr}{$pos}{marker},$stat{$trit}{$chr}{$pos}{mark1},$stat{$trit}{$chr}{$pos}{mark2}),"\n";
		}
	}
}
close OUT;
