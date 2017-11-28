#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
$Script=~s/\.pl//g;
my @Times=localtime();
my $year=$Times[5]+1990;
my $month=$Times[4]+1;
my $day=$Times[3];
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dOut,$fPosi,$Key,$lchr,$mlod);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$mlod,
				"o:s"=>\$dOut,
				"l:s"=>\$lchr,
				"2:s"=>\$fPosi,
				"k:s"=>\$Key,
				) or &USAGE;
&USAGE unless ($mlod and $dOut and $fPosi and $lchr);
open In,$lchr;
my %chrID;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($chr,undef)=split;
	$chrID{$chr}=1;
}
close In;
open In,$fPosi;
my %pos;
my %LG;
my %sca;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,$chr,undef)=split;
	$pos{$id}=$chr;
	if (exists $chrID{$chr}) {
		$LG{$chr}{$id}=1;
	}else{
		$sca{$chr}{$id}=1;
	}
}
close In;
my %PWD;
open In,$mlod;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/MLOD/);
	my ($m1,$m2,$mlod)=split(/\s+/,$_);
	if (!defined $pos{$m1} || !defined $pos{$m2}) {
		die $_;
	}
	next if (!exists $chrID{$pos{$m1}} && !exists $chrID{$pos{$m2}});
	$PWD{$m1}{$m2}=$mlod;
	$PWD{$m2}{$m1}=$mlod;
}
close In;
open Out,">$dOut/$Key.sca.stat";
foreach my $sca (sort keys %sca) {
	my %stat;
	foreach my $marker (sort keys %{$sca{$sca}}) {
		if (!exists $PWD{$marker}) {
			next;
		}
		my @markers=sort {$PWD{$marker}{$b}<=>$PWD{$marker}{$a}} keys %{$PWD{$marker}};
		if (scalar @markers == 0) {
			next;
		}
		$stat{$pos{$markers[0]}}++ if ($PWD{$marker}{$markers[0]} > 5);
	}
	my $target=(sort{$stat{$b} <=> $stat{$a}} keys %stat)[0];
	foreach my $marker (sort keys %{$sca{$sca}}) {
		if (!exists $PWD{$marker}) {
			next;
		}
		$LG{$target}{$marker}=1;
	}
	print Out ">$sca\n";
	foreach my $sca (sort{$stat{$b} <=> $stat{$a}} keys %stat) {
		print Out "$sca:$stat{$sca}\t";
	}
	print Out "\n";
}
foreach my $chr (sort keys %LG) {
	foreach my $marker (sort keys %{$LG{$chr}}) {
		if (!exists $PWD{$marker}) {
			next;
		}
		my @markers=sort{$PWD{$marker}{$b}<=>$PWD{$marker}{$a}} keys %{$PWD{$marker}};
		my $max;
		for (my $i=0;$i<@markers;$i++) {
			if (!exists $PWD{$markers[$i]}) {
				next;
			}
			if ($pos{$markers[$i]} eq $chr) {
				$max=$PWD{$marker}{$markers[$i]};
				last;
			}
		}
		if ($max < 5) {
			print Out "delete:$marker\t$max\n";
			delete $LG{$chr}{$marker};
		}
	}
}
close Out;
open Out,">$dOut/$Key.lg";
foreach my $lg (sort keys %LG) {
	print Out ">$lg\t",scalar keys %{$LG{$lg}},"\n";
	print Out join("\t", keys %{$LG{$lg}}),"\n";
}
close Out;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub USAGE {#
	my $usage=<<"USAGE";
	Program:$Script
	Version:$version	[$month:$day:$year:]
	Contact:Huang Long <huangl\@biomarker.com.cn>
	Options:
		-i	<file>	input mlod file
		-2	<file>	input posi file
		-l	<file>	chr list file
		-o	<dir>	output file
		-k	<str>	output keys of filename
		-h	Help

USAGE
	print $usage;
	exit;
}
