#!/usr/bin/perl 
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
my ($fIn,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my @xls=glob("$fIn/*.xls");
my %stat;
foreach my $xls (@xls) {
	open In,$xls;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($laneID,$libID,$project,$sampleID,$majorID,$rawreads,$cleanreads,$rawbase,$cleanbase,$rawQ30,$cleanQ30,undef,undef,$rawGC,$cleanGC,$rawdup,$cleandup,$adapter,$pollution1,,$pollution2,undef)=split(/\s+/,$_);
		$stat{$project}{$sampleID}{rawreads}+=$rawreads;
		$stat{$project}{$sampleID}{cleanreads}+=$cleanreads;
		$stat{$project}{$sampleID}{rawbase}+=$rawbase;
		$stat{$project}{$sampleID}{cleanbase}+=$cleanbase;
		$stat{$project}{$sampleID}{rawQ30}+=$rawbase*$rawQ30;
		$stat{$project}{$sampleID}{cleanQ30}+=$cleanbase*$cleanQ30;
		$stat{$project}{$sampleID}{rawGC}+=$rawbase*$rawGC;
		$stat{$project}{$sampleID}{cleanGC}+=$cleanbase*$cleanGC;
		$stat{$project}{$sampleID}{rawdup}+=$rawbase*$rawGC;
		$stat{$project}{$sampleID}{cleandup}+=$cleanbase*$cleanGC;
		$stat{$project}{$sampleID}{adapter}+=$adapter*$cleanGC;
	}
	close In;
}
mkdir $fOut if (!-d $fOut);
foreach my $project (sort keys %stat) {
	open Report,">$fOut/$project.report.xls";
	print Report join("\t","#sampleID","RawBase","CleanBase","Q30(%)","GC(%)"),"\n";
	foreach my $sampleID (sort keys %{$stat{$project}}) {
		my @out;
		push @out,$sampleID;
		push @out,$stat{$project}{$sampleID}{rawbase};
		push @out,$stat{$project}{$sampleID}{cleanbase};
		push @out,sprintf("%.2f",$stat{$project}{$sampleID}{cleanQ30}/$stat{$project}{$sampleID}{cleanbase});
		push @out,sprintf("%.2f",$stat{$project}{$sampleID}{cleanGC}/$stat{$project}{$sampleID}{cleanbase});
		print Report join("\t",@out),"\n";
	}
	close Report;
}
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
	-i	<file>	input file 
	-o	<file>	output file
USAGE
	print $usage;
	exit;
}
