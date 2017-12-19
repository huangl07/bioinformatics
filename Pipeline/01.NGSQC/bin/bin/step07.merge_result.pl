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
open In,"$dIn/06.QCreport/10.report.xls";
my %stat;
my ($project,$laneID,$libID);
while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($laneID,$libID,$project,$majorID,$sampleID,$rawdata,$cleandata,$cleanreads,$cleanper,$cleanQ30,$cleanGC,$cleandup,$pollution1,$pollution2,$pullution3)=split(/\t/,$_);
#		$stat{$project}{$sampleID}{rawdata}+=$rawdata;
#		$stat{$project}{$sampleID}{cleandata}+=$cleandata;
#		$stat{$project}{$sampleID}{cleanreads}+=$cleanreads;
#		$stat{$project}{$sampleID}{cleanQ30}+=$cleandata*$cleanQ30;
#		$stat{$project}{$sampleID}{cleanGC}+=$cleandata*$cleanGC;
#		$stat{$project}{$sampleID}{cleandup}+=$cleanbase*$cleanGC;
#		$stat{$project}{$sampleID}{adapter}+=$adapter*$cleanGC;
		#$stat{$project}{$sampleID}{$laneID}{pollution1}+=$pollution1;
		#$stat{$project}{$sampleID}{$laneID}{pollution2}+=$pollution2;
		$stat{$project}{$laneID}{$libID}{$sampleID}{rawdata}=$rawdata;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleandata}=$cleandata;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleanreads}=$cleanreads;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleanper}=$cleanper;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleanQ30}=$cleanQ30;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleanGC}=$cleanGC;
		$stat{$project}{$laneID}{$libID}{$sampleID}{cleandup}=$cleandup;
		$stat{$project}{$laneID}{$libID}{$sampleID}{pollution1}=$pollution1;
		$stat{$project}{$laneID}{$libID}{$sampleID}{pollution2}=$pollution2;
		$stat{$project}{$laneID}{$libID}{$sampleID}{pullution3}=$pullution3;
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
	open FQ,">$dOut/$project/fq-$project.list";
#	print FQ join("\t",
	open Report,">$dOut/$project/$project.report.xls";
	print Report join("\t","#sampleID","RawBase","CleanBase","CleanReads","CleanBasePer(%)","Q30(%)","GC(%)","Pollution-species","Pollution-rate(%)"),"\n";
	foreach my $sampleID (sort keys %{$stat{$project}}) {
		open FQ,">$dOut/$project/fq-$project.list";
		print FQ join("\t",$sampleID,"$dIn/03.CleanData/*$project*$sampleID*clean.1.fastq.gz","$dIn/03.CleanData/*$project*$sampleID*clean.2.fastq.gz"),"\n"; 
		my @out;
		push @out,$sampleID;
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{rawdata};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleandata};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleanreads};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleanper};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleanQ30};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleanGC};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{cleandup};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{pollution1};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{pollution2};
		push @out,$stat{$project}{$laneID}{$libID}{$sampleID}{pullution3};
		print Report join("\t",@out),"\n";
	}
	
	close Report;
	close FQ;
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
