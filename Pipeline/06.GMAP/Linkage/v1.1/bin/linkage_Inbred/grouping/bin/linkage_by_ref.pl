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
my ($fIn,$fOut,$fPosi,$Key);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"2:s"=>\$fPosi,
				"k:s"=>\$Key,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $fPosi and $Key);
mkdir $fOut if (!-d $fOut);
open In,$fIn;
my %Marker;
my $Head=<In>;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /ID/);
	my @temp=split(/\s+/,$_);
	$Marker{$temp[0]}=$_;
}
close In;
open In,$fPosi;
my %Posi;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ );
	my @temp=split(/\s+/,$_);
	next if (!exists $Marker{$temp[0]});
	my $Gm=$temp[1];
	$Posi{$Gm}{$temp[0]}=$temp[2];
}
close In;
open Out,">$fOut/$Key.lg";
foreach my $LG (sort keys %Posi) {
	print Out ">$LG\t",scalar keys %{$Posi{$LG}},"\n",join("\t",keys %{$Posi{$LG}}),"\n";
	open Sort,">$fOut/$Key.$LG.sort";
	my @Marker=sort {$Posi{$LG}{$a}<=>$Posi{$LG}{$b}} keys %{$Posi{$LG}};
	my $LG1 = $LG;
	$LG1=~s/\D+//g;
	print Sort "group $LG1\n";
	foreach my $pos (@Marker) {
		print  Sort $pos,"\t",$Posi{$LG}{$pos},"\n";
	}
	close Sort;
	open Geno,">$fOut/$Key.$LG.genotype";
	print Geno $Head;
	foreach my $m (@Marker) {
		print Geno $Marker{$m},"\n";
	}
	close Geno;
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
		-i	<file>	input genotype file
		-2	<file>	input posi file
		-o	<dir>	output file
		-k	<str>	output keys of filename
		-h	Help

USAGE
	print $usage;
	exit;
}
