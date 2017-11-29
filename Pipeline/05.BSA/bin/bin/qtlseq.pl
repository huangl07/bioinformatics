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
my ($fIn,$fOut,$PID,$BID,$Pdep,$Bdep,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$fIn,
				"out:s"=>\$fOut,
				"pid:s"=>\$PID,
				"bid:s"=>\$BID,
				"pdep:s"=>\$Pdep,
				"Bdep:s"=>\$Bdep,
				"popt:s"=>\$popt,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $BID);
$Pdep||=10;
$Bdep||=10;
$PID||="";
$popt||="F2";
my %Indi;
if ($PID ne "") {
	my ($P1,$P2)=split(/\,/,$PID);
	$P2||=$P1;
	$Indi{$P1}="P1";
	$Indi{$P2}="P2";
}
my ($B1,$B2)=split(/\,/,$BID);
$Indi{$B1}="B1";
$Indi{$B2}="B2";
if ($fIn =~ /gz/) {
	open In,"gunzip -dc $fIn|";
}else{
	open In,$fIn;
}
open Out,">$fOut";
print Out join("\t","#chr","pos","ref","alt","anno","P1GT","P1AD","P2GT","P2AD","B1GT","B1AD","B2GT","B2AD","INDEX1","INDEX2","DELTA"),"\n";
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
		my $ann="$ids";
		if($info=~/ANN=([^\;]*)/){$ann=$1;my @ann=split(/\|/,$ann);$ann=join(/|/,$ann[1],$ann[2],$ann[3],$ann[4])}
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
		next if ($info{B1}{gt}  eq "./." || $info{B2}{gt} eq "./." ||$info{B1}{dp} < $Bdep || $info{B2}{dp} < $Bdep);
		my @b1=split(/\/|\|/,$info{B1}{gt});
		my @b2=split(/\/|\|/,$info{B2}{gt});
		my @ad1=split(/\,/,$info{B1}{ad});
		my @ad2=split(/\,/,$info{B2}{ad});
		my %stat;
		$stat{$b1[0]}++;
		$stat{$b1[1]}++;
		$stat{$b2[0]}++;
		$stat{$b2[1]}++;
		my @geno=sort {$stat{$a}<=>$stat{$b}} keys %stat;
		my $index1;
		my $index2;
		my $delta;
		if ($popt eq "F2") {
			next if (scalar @geno!=2);
			if ($PID eq "") {
				$index1=$ad1[$geno[0]]/$info{B1}{dp};
				$index2=$ad2[$geno[0]]/$info{B2}{dp};
				$delta=abs($index1-$index2);
			}else{
				my ($p1,$p2)=split(/\/|\|/,$info{P1}{gt});
				if (!exists $info{P2}) {
					if ($p1 eq $geno[0]) {
						$info{P2}{gt}="$geno[0]q/$geno[0]";
						$info{P2}{dp}=10;
					}else{
						$info{P2}{gt}="$geno[1]/$geno[1]";
						$info{P2}{dp}=10;
					}
				}
				my ($p3,$p4)=split(/\/|\|/,$info{P2}{gt});
				next if ($p1 ne $p2 || $p3 ne $p4);
				next if	($info{P2}{dp}< $Pdep || $info{P1}{dp} < $Bdep);
				next if ($p1 eq $p3);
				next if (!exists $stat{$p1} || !exists $stat{$p3});
				$index1=$ad1[$p1]/$info{B1}{dp};
				$index2=$ad2[$p1]/$info{B2}{dp};
				$delta=$index1-$index2;
			}
		}else{
			next if (scalar @geno!=2);
			my ($p1,$p2)=split(/\/|\|/,$info{P1}{gt});
			my ($p3,$p4)=split(/\/|\|/,$info{P2}{gt});
			next if ($p1 eq $p2 || $p3 ne $p4);
			next if ($p1 ne $p2 || $p3 eq $p4);
			next if	($info{P2}{dp}< $Pdep || $info{P1}{dp} < $Bdep);
			if ($p3 eq $p1 && $p1 eq $p2) {#nnxnp
				$index1=$ad1[$p1]/$info{B1}{dp};
				$index1=$ad2[$p1]/$info{B2}{dp};
				$delta=$index1-$index2;
			}elsif ($p3 eq $4 && $p3 eq $p1) {#lmxll
				$index1=$ad1[$p2]/$info{B1}{dp};
				$index2=$ad2[$p2]/$info{B1}{dp};
				$delta=$index1-$index2;
			}else{
				next;
			}
		}
		$info{B1}{gt}||="--";
		$info{B1}{ad}||="--";
		$info{B2}{gt}||="--";
		$info{B2}{ad}||="--";
		$info{P1}{gt}||="--";
		$info{P1}{ad}||="--";
		$info{P2}{gt}||="--";
		$info{P2}{ad}||="--";
		print Out join("\t",$chr,$pos,$ref,$alt,$ann,$info{P1}{gt},$info{P1}{ad},$info{P2}{gt},$info{P2}{ad},$info{B1}{gt},$info{B2}{ad},$info{B1}{gt},$info{B2}{ad},$index1,$index2,$delta),"\n";
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
	-popt	<str>	population type default F2
USAGE
	print $usage;
	exit;
}
