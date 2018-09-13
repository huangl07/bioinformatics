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
my ($dIn,$list,$dOut);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$dIn,
				"s:s"=>\$list,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($dIn and $list and $dOut);
open LIST,"$list";
my @id;
while (<LIST>){
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	push @id,$_;
}
close LIST;

open In,"$dIn/06.QCreport/qc-report.xls";
#open In,"$dIn/06.QCreport/6.report.xls";
my %stat;
my ($project,$laneID,$libID);
while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($laneID,$libID,$project,$majorID,$sampleID,$rawread,$cleanread,$rawdata,$cleandata,$rawq30,$cleanq30,$rawgc,$cleangc,$rawdup,$cleandup,$ada,$pollution1,$pollution2,$pollution3)=split(/\s+/,$_);
		$stat{$project}{$sampleID}{majorID}=$majorID;
		$stat{$project}{$sampleID}{rawdata}+=$rawdata;
		$stat{$project}{$sampleID}{cleandata}+=$cleandata;
		$stat{$project}{$sampleID}{cleanreads}+=$cleanread;
		$stat{$project}{$sampleID}{cleanQ30}+=$cleandata*$cleanq30;
		$stat{$project}{$sampleID}{cleanGC}+=$cleandata*$cleangc;
		$stat{$project}{$sampleID}{dup}+=$cleandata*$cleandup;
		$stat{$project}{$sampleID}{pollution1}=$pollution1;
		$stat{$project}{$sampleID}{pollution2}=$pollution2;
		$stat{$project}{$sampleID}{pollution3}=$pollution3;
}
close In;
open Out,">$dIn/06.QCreport/6.report.xls";
print Out join("\t","#ProjectID\tMajorID\tSampleID","Raw Data(bp)","Clean Data(bp)","Clean Reads(reads)","CleanBasePer(%)","Q30(%)","GC(%)","Dup(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my$id(@id){
	foreach my $project (sort keys %stat){
		foreach my $sampleID (sort keys %{$stat{$project}}) {
			if($stat{$project}{$sampleID}{majorID} eq $id){
				print Out join("\t",$project,$stat{$project}{$sampleID}{majorID},$sampleID,$stat{$project}{$sampleID}{rawdata},$stat{$project}{$sampleID}{cleandata},$stat{$project}{$sampleID}{cleanreads},sprintf("%.2f",100*$stat{$project}{$sampleID}{cleandata}/$stat{$project}{$sampleID}{rawdata}),sprintf("%.2f",$stat{$project}{$sampleID}{cleanQ30}/$stat{$project}{$sampleID}{cleandata}),sprintf("%.2f",$stat{$project}{$sampleID}{cleanGC}/$stat{$project}{$sampleID}{cleandata}),sprintf("%.2f",$stat{$project}{$sampleID}{dup}/$stat{$project}{$sampleID}{cleandata}),$stat{$project}{$sampleID}{pollution1},$stat{$project}{$sampleID}{pollution2},$stat{$project}{$sampleID}{pollution3}),"\n";
			}
		}
	}
}
close Out;
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
	`grep "$project" $dIn/03.CleanData/fastq.list >$dOut/$project/fq.list`;
	open Report,">$dOut/$project/$project.report.xls";
	print Report join("\t","#SampleID","CleanReads(reads)","CleanBase(bp)","Q30(%)","GC(%)","Pollution-species","Pollution-rate(%)"),"\n";
	open In,"$dOut/$project/fq.list";
	open Out,">$dOut/$project/fq-$project.list";
	while (<In>) {
		chomp;
		next if ($_ eq ""|| /^$/);
		my($sample,$fq1,$fq2,$type)=split/\t/;
		my $name=(split(/\:/,$sample))[-1];
		print Out "$name\t$fq1\t$fq2\t$type\n";
	}
	close In;
	close Out;
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
	-s	<file>	input sample.list
	-o	<dir>	output dir
USAGE
	print $usage;
	exit;
}
