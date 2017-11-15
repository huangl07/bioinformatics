#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$report,$fmtrix);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"m:s"=>\$fmtrix,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $fmtrix);
my %TsTv=(
	"AT"=>"Tv",
	"AG"=>"Ts",
	"AC"=>"Tv",
	"GT"=>"Tv",
	"CT"=>"Ts",
	"CG"=>"Tv",
);
open In,$fIn;
my @indi;
my %vcfstat;
my %diff;
my %pop;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	next if ($Filter ne "PASS"  && $Filter ne "FILTER");
	my @Ale=split(/\,/,join(",",$ref,$alt));
	if (scalar @indi ==0) {
		push @indi,@geno;
		next;
	}else{
		my @Geno;
		my %geno;
		my %nomiss;
		#$stat{$ref}++;
		for (my $i=0;$i<@geno;$i++) {
			my ($gt,$ad,$dp,undef)=split(/\:/,$geno[$i]);
			push @Geno,$gt;
			$geno{$gt}++;
			my ($g1,$g2)=split(/\//,$gt);
			if ($gt eq "./." || $Ale[$g1] eq "*" || $Ale[$g2] eq "*"){
				$vcfstat{$indi[$i]}{miss}++;
				next;
			}
			$nomiss{$g1}++;
			if ($gt eq "0/0"){
				$vcfstat{$indi[$i]}{ref}++;
				next;
			}
			$vcfstat{$indi[$i]}{total}++;
			$vcfstat{$indi[$i]}{homo}++ if($g1 eq $g2);
			my %change;
			$change{$Ale[$g1]}=1;
			$change{$Ale[$g2]}=1;
			$change{$ref}=1;
			if (scalar keys %change >= 3) {
				$vcfstat{$indi[$i]}{multi}++;
			}else{
				my $change=join("",sort keys %change);
				$vcfstat{$indi[$i]}{$TsTv{$change}}++;
			}
		}
		if (scalar keys %nomiss > 1) {
			$pop{total}++;
			my %stat;
			foreach my $gt (keys %geno) {
				my ($g1,$g2)=split(/\//,$gt);
				if ($gt eq "./." || $Ale[$g1] eq "*" || $Ale[$g2] eq "*"){
					$pop{miss}++;
					next;
				}
				if ($gt eq "0/0"){
					$pop{ref}++;
				}
				$pop{homo}++ if($g1 eq $g2);
				$stat{$g1}++;
				$stat{$g2}++;
			}
			my ($g1,$g2,undef)=sort{$stat{$b}<=>$stat{$a}} keys %stat;
			my $change=join("",sort ($Ale[$g1],$Ale[$g2]));
			$pop{$TsTv{$change}}++;
		}

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
open Out,">$fmtrix";
my @sample=sort keys %vcfstat;
print Out join("\t","",@sample),"\n";
for (my $i=0;$i<@sample;$i++) {
	my @out;push @out,$sample[$i];
	for (my $j=0;$j<@sample;$j++) {
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
print Out "#sampleID\tSNPnumber\tTransition\tTransvertion\tTs/Tv\tHeterozygosity Number\tHomozygosity Number\tMiss Number\tRef Number\n";
foreach my $sample (sort keys %vcfstat) {
	print Out join("\t",$sample,$vcfstat{$sample}{total},$vcfstat{$sample}{Ts},$vcfstat{$sample}{Tv},sprintf("%.2f",$vcfstat{$sample}{Ts}/$vcfstat{$sample}{Tv}),$vcfstat{$sample}{total}-$vcfstat{$sample}{homo},$vcfstat{$sample}{homo},$vcfstat{$sample}{miss},$vcfstat{$sample}{ref}),"\n";
}
$pop{miss}||=0;
$pop{Ts}||=0;
$pop{ref}||=0;
print Out join("\t","pop",$pop{total},$pop{Ts},$pop{Tv},sprintf("%.2f",$pop{Ts}/$pop{Tv}),$pop{total}*(scalar @sample)-$pop{homo},$pop{homo},$pop{miss},$pop{ref}),"\n";
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
  -o	<file>	output file name
  -m	<file>	output metric file
  -h         Help

USAGE
        print $usage;
        exit;
}
