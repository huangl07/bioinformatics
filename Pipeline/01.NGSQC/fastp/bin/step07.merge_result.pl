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
my ($dIn,$dOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$dIn,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($dIn and $dOut);
open In,"$dIn/06.QCreport/6.report.xls";
my %stat;
my ($project,$laneID,$libID);
while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($laneID,$libID,$project,$majorID,$sampleID,$rawdata,$cleandata,$cleanreads,$cleanper,$cleanQ30,$cleanGC,$dup,$pollution1,$pollution2,$pollution3)=split(/\s+/,$_);
#		$stat{$project}{$sampleID}{rawdata}+=$rawdata;
		$stat{$project}{$sampleID}{cleandata}+=$cleandata;
		$stat{$project}{$sampleID}{cleanreads}+=$cleanreads;
		$stat{$project}{$sampleID}{cleanQ30}+=$cleandata*$cleanQ30;
		$stat{$project}{$sampleID}{cleanGC}+=$cleandata*$cleanGC;
		$stat{$project}{$sampleID}{pollution1}=$pollution1;
		$stat{$project}{$sampleID}{pollution2}=$pollution2;
		$stat{$project}{$sampleID}{pollution3}=$pollution3;
}
close In;


mkdir "$dOut/" if (!-d "$dOut/");
#mkdir "$dOut/$project/Figure" if (!-d "$dOut/$project/Figure");
#`ln -s $dIn/fig/*$project* $dOut/$project/Figure`;

foreach  $project (sort keys %stat) {
	`mkdir "$dOut/$project" ` if (!-d "$dOut/$project");
	`mkdir "$dOut/$project/Figure"` if (!-d "$dOut/$project/Figure");
	`ln -s $dIn/06.QCreport/fig/*$project* $dOut/$project/Figure`;
	`mkdir "$dOut/$project/cleandata" ` if (!-d "$dOut/$project/cleandata");
	`ln -s $dIn/03.CleanData/*$project*clean*fastq.gz  $dOut/$project/cleandata`;
	`mkdir "$dOut/$project/rawdata" ` if (!-d "$dOut/$project/rawdata");
	`ln -s $dIn/01.RawData/*$project* $dOut/$project/rawdata`;
	`grep "$project" $dIn/03.CleanData/fastq.list >$dOut/$project/fq-$project.list`;
	open Report,">$dOut/$project/$project.report.xls";
	print Report join("\t","#SampleID","CleanReads(reads)","CleanBase(bp)","Q30(%)","GC(%)","Pollution-species","Pollution-rate(%)"),"\n";
	foreach my $sampleID (sort keys %{$stat{$project}}) {
		print Report join("\t",$sampleID,$stat{$project}{$sampleID}{cleanreads},$stat{$project}{$sampleID}{cleandata},sprintf("%.2f",$stat{$project}{$sampleID}{cleanQ30}/$stat{$project}{$sampleID}{cleandata}),sprintf("%.2f",$stat{$project}{$sampleID}{cleanGC}/$stat{$project}{$sampleID}{cleandata}),$stat{$project}{$sampleID}{pollution1},$stat{$project}{$sampleID}{pollution2},$stat{$project}{$sampleID}{pollution3}),"\n";
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
Contact: caixia.tian

Usage:
  Options:
	-i	<dir>	input dir
	-o	<dir>	output dir
USAGE
	print $usage;
	exit;
}
