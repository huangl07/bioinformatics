#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$PID,$BID,$Pdep,$Bdep);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$fIn,
				"out:s"=>\$fOut,
				"pid:s"=>\$PID,
				"bid:s"=>\$BID,
				"pdep:s"=>\$Pdep,
				"Bdep:s"=>\$Bdep,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $BID);
$Pdep||=10;
$Bdep||=10;
$PID||="";
my %Indi;
if ($PID ne "") {
	my ($P1,$P2)=split(/\,/,$PID);
	$P2||=$P1;
	$Indi{$P1}="P1";
	if ($P1 ne $P2) {
		$Indi{$P2}="P2";
	}
}
$Indi{$BID}="B";
if ($fIn =~ /gz/) {
	open In,"gunzip -dc $fIn|";
}else{
	open In,$fIn;
}
open Out,">$fOut";
print Out join("\t","#chr","pos","ref","alt","anno","PID","PAD","BID","BAD","INDEX"),"\n";
my @indi;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	my ($chr,$pos,$ids,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\s+/,$_);
	if (/^#/) {
		push @indi,@geno;
	}else{
		my %info;
		my @format=split(/:/,$format);
		my $ann=$ids;
		if($info=~/ANN=([^\;]*)/){$ann=$1;my @ann=split(/\|/,$ann);$ann=join("|",$ann[1],$ann[2],$ann[3],$ann[4])}
		for (my $i=0;$i<@indi;$i++) {
			next if (!exists $Indi{$indi[$i]});
			my $id=$Indi{$indi[$i]};
			my @info=split(/:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				$info{$id}{gt}=$info[$j] if ($format[$j] eq "GT");
				$info{$id}{ad}=$info[$j] if ($format[$j] eq "AD");
				$info{$id}{dp}=$info[$j] if ($format[$j] eq "DP");
			}
		}
		if (!exists $info{P1}) {
			$info{P1}{gt}="0/0";
			$info{P1}{dp}="10";
			$info{P1}{ad}="10,0";
		}
		next if ($info{P1}{gt} eq "./." || $info{B}{gt} eq "./.");
		next if ($info{P1}{gt} eq $info{B}{gt});
		next if ($info{P1}{dp} < $Pdep || $info{B}{dp} < $Bdep);
		my ($p1,$p2)=split(/\/|\|/,$info{P1}{gt});
		my ($b1,$b2)=split(/\/|\|/,$info{B}{gt});
		next if ($p1 ne $p2);
		next if ($p1 ne $b1 && $p1 ne $b2 && $b1 ne $b2);
		my $mut=$b1;
		if ($p1 == $b1) {
			$mut=$b2;
		}
		if (exists $info{P2} && exists $info{P1} ) {
			my ($p3,$p4)=split(/\/|\|/,$info{P2}{gt});
			next if ($p3 ne $p4);
			next if ($info{P1}{gt} ne $info{P2}{gt});
			next if($p3 ne $b1 && $p3 ne $b2);
			$mut=$p3;
		}
		my @dp=split(/\,/,$info{B}{ad});
		my $sum=0;
		if ($b1 eq $b2) {
			$sum=$dp[$b1];
		}else{
			$sum=$dp[$b1]+$dp[$b2];
		}
		my $snpindex=$dp[$mut]/$sum;
		$info{P1}{gt}||="--";
		$info{P1}{ad}||="--";
		$info{B}{gt}||="--";
		$info{B}{ad}||="--";
		print Out join("\t",$chr,$pos,$ref,$alt,$ann,$info{P1}{gt},$info{P1}{ad},$info{B}{gt},$info{B}{ad},$snpindex),"\n";
	}
}
close In;
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-vcf	<file>	input file 
	-out	<file>	output file
	-pid	<str>	wild parent id may not
	-bid	<str>	mut bulk id
	-pdep	<str>	parent depth default 10
	-bdep	<str>	mut build depth default 10
USAGE
	print $usage;
	exit;
}
