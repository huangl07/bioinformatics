#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fPos,$paternalID,$maternalID,$ParentDepth,$OffDepth);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"p:s"=>\$fPos,
	"PID:s"=>\$paternalID,
	"MID:s"=>\$maternalID,
	"Pdep:s"=>\$ParentDepth,
	"Odep:s"=>\$OffDepth,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $paternalID and $maternalID);
$ParentDepth||=10;
$OffDepth||=10;
open In,$fIn;
open Out,">$fOut";
open Pos,">$fPos";
my @indi;
my $PID;
my $MID;
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
		my @outid;
		for (my $i=0;$i<@indi;$i++) {
			if ($indi[$i] eq $paternalID) {
				$PID=$i;
				next;
			}
			if ($indi[$i] eq $maternalID) {
				$MID=$i;
				next;
			}
			push @outid,$indi[$i];
		}
		die "error input PID or MID" if (!defined $PID || !defined $MID);
		print Out join("\t","#MarkerID","type",@outid),"\n";
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
		$n++;
		my @out;
		push @out,"Marker$n";
		my ($PAG,$PAD,$PTD,undef)=split(/\:/,$geno[$PID]);
		my ($MAG,$MAD,$MTD,undef)=split(/\:/,$geno[$MID]);
		next if ($geno[$PID] eq "./." || $geno[$MID] eq "./.");
		my ($p1,$p2)=split(/\//,$PAG);
		my ($m1,$m2)=split(/\//,$MAG);
		my ($pd1,$pd2)=split(/\,/,$PAD);
		my ($md1,$md2)=split(/\,/,$MAD);
		next if ($MAG eq $PAG && $p1 eq $p2);
		next if ($PAG eq "./." || $MAG eq "./.");
		next if ($pd1+$pd2 < $ParentDepth || $md1 + $md2 < $ParentDepth);
		my %geno;
		$geno{$p1}++;
		$geno{$p2}++;
		$geno{$m1}++;
		$geno{$m2}++;
		my %ale;
		if (scalar keys %geno == 4) { #abxcd
			 %ale=(
				$p1=>"a",
				$p2=>"b",
				$m1=>"c",
				$m2=>"d",
			);
			push @out,"abxcd";
		}elsif (scalar keys %geno == 3 && $p1 eq $p2 && $m1 ne $m2) {#ccxab
			 %ale=(
				$p1=>"c",
				$m1=>"a",
				$m2=>"b",
			);
			push @out,"ccxab";
		}elsif (scalar keys %geno == 3 && $p1 ne $p2 && $m1 eq $m2) {#abxcc
			 %ale=(
				$p1=>"a",
				$p2=>"b",
				$m1=>"c",
			);
			push @out,"abxcc";
		}elsif (scalar keys %geno == 3 && $p1 ne $p2 && $m1 ne $m2) {#efxeg
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"f",
				$ale[1]=>"g",
				$ale[2]=>"e",
			);
			push @out,"efxeg";
		}elsif (scalar keys %geno == 2 && $PAG eq $MAG) {#hkxhk
			%ale=(
				$p1=>"h",
				$p2=>"k",
			);
			push @out,"hkxhk";
		}elsif (scalar keys %geno == 2 && $p1 eq $p2 && $m1 ne $m2) {#nnxnp
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"p",
				$ale[1]=>"n",
			);
			push @out,"nnxnp";
		}elsif (scalar keys %geno == 2 && $p1 ne $p2 && $m1 eq $m2) {#lmxll
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"m",
				$ale[1]=>"l",
			);
			push @out,"lmxll";
		}elsif (scalar keys %geno == 2 && $p1 eq $p2 && $m1 eq $m2) {#aaxbb
			 %ale=(
				$p1=>"a",
				$m1=>"b",
			);
			push @out,"aaxbb";
		}
		for (my $i=0;$i<@indi;$i++) {
			next if ($i eq $PID || $i eq $MID);
			my ($GAG,$GAD,$GTD,undef)=split(/\:/,$geno[$i]);
			if ($GAG eq "./.") {
				push @out,"--";
				next;
			}
			my ($od1,$od2)=split(/\,/,$GAD);
			if ($od1+$od2 < $OffDepth) {
				push @out,"--";
			}else{
				my ($g1,$g2)=split(/\//,$GAG);
				if (!exists $ale{$g1} || !exists $ale{$g2}) {
					push @out,"--";
				}else{
					my $geno=join("",sort ($ale{$g1},$ale{$g2}));
					push @out,"$geno";
				}
			}
		}
		print Out join("\t",@out),"\n";
		print Pos join("\t",$out[0],$chr,$pos),"\n";
	}
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT,@INDI)=split;
}
close Out;
close In;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -i	<file>	input file vcf file
  -o	<file>	output marker file name
  -p	<file>	output pos file name
  -PID	<str>	paternal ID
  -MID	<str>	maternale ID
  -Pdep	<num>	parentDepth to filter
  -Odep	<num>	offspringDepth to filter
  -h         Help

USAGE
        print $usage;
        exit;
}
