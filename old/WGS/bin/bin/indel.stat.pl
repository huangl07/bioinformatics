#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$report,$fmtrix,$flen);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"m:s"=>\$fmtrix,
	"l:s"=>\$flen,
			) or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
open Flen,">$flen";
my @indi;
my %vcfstat;
my %diff;
my %len;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	next if ($Filter ne "PASS" && $Filter ne "SNP" && $Filter ne "INDEL" && $Filter ne "FILTER");
	my @Ale=split(/\,/,join(",",$ref,$alt));
	if (scalar @indi ==0) {
		push @indi,@geno;
		next;
	}else{
		my @Geno;
		my %geno;
		my %pchange;
		my $sum=0;
		my $sumlen=0;
		for (my $i=0;$i<@geno;$i++) {
			my ($gt,$ad,$dp,undef)=split(/\:/,$geno[$i]);
			push @Geno,$gt;
			my ($g1,$g2)=split(/\//,$gt);
			if ($gt eq "./." || $Ale[$g1] eq "*" || $Ale[$g2] eq "*"){
				$vcfstat{$indi[$i]}{miss}++;
				next;
			}
			if ($gt eq "0/0"){
				$vcfstat{$indi[$i]}{ref}++;
				next;
			}
			$geno{$gt}++ if($gt ne "./.");
			$sum+=$dp if($gt ne "./.");
			if (length($ref)*2 < length($Ale[$g1])+length($Ale[$g2])) {
				$vcfstat{$indi[$i]}{delete}++;
			}else{
				$vcfstat{$indi[$i]}{insert}++;
			}
			$vcfstat{$indi[$i]}{total}++;
			$vcfstat{$indi[$i]}{dp}+=$dp;
			$vcfstat{$indi[$i]}{homo}++ if($g1 eq $g2);
			$vcfstat{$indi[$i]}{lenth}+=(length($Ale[$g1])+length($Ale[$g2]))/2;
			my $len=length($Ale[0]);
			my $len1=length($Ale[$g1])-length($Ale[0]);#-length($Ale[0]);
			my $len2=length($Ale[$g2])-length($Ale[0]);
			$sumlen=$len*2-$len1-$len2;
			print Flen $indi[$i],"\t",$len1-$len,"\n" if ($len1 != $len);
			print Flen $indi[$i],"\t",$len1-$len,"\n" if ($len2 != $len);
		}
		$vcfstat{pop}{dp}+=$sum if(scalar keys %geno != 1);
		$vcfstat{pop}{len}+=$sumlen if(scalar keys %geno != 1);
		 if(scalar keys %geno == 2){
			if ($sumlen < 0) {
				$vcfstat{pop}{insert}++;
			}else{
				$vcfstat{pop}{delete}++;
			}
		 };
		 $vcfstat{pop}{total}++ if(scalar keys %geno !=1);
		for (my $i=0;$i<@Geno;$i++) {
			for (my $j=$i+1;$j<@Geno;$j++) {
				if($Geno[$i] ne $Geno[$j] && $Geno[$i] ne "./." && $Geno[$j] ne "./."){
					$diff{$indi[$j]}{$indi[$i]}++ ;
					$diff{$indi[$i]}{$indi[$j]}++ ;
				};
			}
		}
	}
}
close In;
close Flen;
open Out,">$fmtrix";
my @sample=sort keys %vcfstat;
print Out join("\t","",@sample),"\n";
for (my $i=0;$i<@sample;$i++) {
	next if ($sample[$i] eq "pop");
	my @out;push @out,$sample[$i];
	for (my $j=0;$j<@sample;$j++) {
		next if ($sample[$j] eq "pop");
		if (!exists $diff{$sample[$i]}{$sample[$j]}) {
			push @out,0;
		}else{
			push @out,$diff{$sample[$i]}{$sample[$j]};
		}
	}
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut";
print Out "Sample ID\tInsert Number\tDelete Number\tHeterozygosity Number\tHomozygosity Number\tAverage Length\tAverage Depth\tMiss Number\tRef Number\n";
foreach my $sample (sort keys %vcfstat) {
	$vcfstat{$sample}{lenth}||=0;
	next if($sample eq "pop");
	print Out join("\t",$sample,$vcfstat{$sample}{insert},$vcfstat{$sample}{delete},$vcfstat{$sample}{homo},$vcfstat{$sample}{total}-$vcfstat{$sample}{homo},sprintf("%.2f",$vcfstat{$sample}{lenth}/$vcfstat{$sample}{total}),sprintf("%.2f",$vcfstat{$sample}{dp}/$vcfstat{$sample}{total}),$vcfstat{$sample}{miss},$vcfstat{$sample}{ref}),"\n";
}
print Out join("\t","pop",$vcfstat{pop}{insert},$vcfstat{pop}{delete},"--","--",sprintf("%.2f",$vcfstat{pop}{lenth}/$vcfstat{pop}{total}),sprintf("%.2f",$vcfstat{pop}{dp}/$vcfstat{pop}{total}),"--","--"),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output stat name
  -m	<file>	output matrix file
  -l	<file>	output len file
  -h         Help

USAGE
        print $usage;
        exit;
}
