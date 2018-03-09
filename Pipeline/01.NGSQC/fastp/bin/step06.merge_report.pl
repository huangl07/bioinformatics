#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($rawdup , $cleandup , $qc , $dOut,$pollution);
GetOptions(
				"help|?" =>\&USAGE,
				"rawdup:s"=>\$rawdup,
				"cleandup:s"=>\$cleandup,
				"pollution:s"=>\$pollution,
				"qc:s"=>\$qc,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($rawdup and $qc and $cleandup and $pollution and $dOut);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
my %stat;
open In,"$qc/stat.list";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$stats)=split(/\s+/,$_);
	open Stat,$stats;
	while (<Stat>) {
		next if ($_ eq "" ||/^$/||/#/);
		($sample,$stat{$sample}{rawread},$stat{$sample}{rawbase},$stat{$sample}{rawq20},$stat{$sample}{rawq30},$stat{$sample}{rawGC},$stat{$sample}{ada},$stat{$sample}{cleanread},$stat{$sample}{cleanbase},$stat{$sample}{cleanq20},$stat{$sample}{cleanq30},$stat{$sample}{cleanGC})=split(/\t/,$_);
	}
	close Stat;
}
close In;
open In,$rawdup;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$dup)=split(/\s+/,$_);
	open Dup,$dup;
	while (<Dup>) {
		chomp;
		next if ($_ eq "" || /^$/);
		if (/Percentage/) {
			my $rdup=(split(/\s+/,$_))[2];
			$stat{$sample}{rdup}=$rdup;
		}
	}
	close Dup;
}
close In;
open In,$cleandup;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$dup)=split(/\s+/,$_);
	open Dup,$dup;
	while (<Dup>) {
		chomp;
		next if ($_ eq "" || /^$/);
		if (/Percentage/) {
			$stat{$sample}{cdup}=(split(/\s+/,$_))[2];
		}
	}
	close Dup;
}
close In;
open In,$pollution;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$po1,$po2,$po3)=split(/\s+/,$_);
	$stat{$sample}{pollutions}=$po1;
	$stat{$sample}{pollutionr}=$po2;
	$stat{$sample}{pollutionn}=$po3;
}
close In;
open Out,">$dOut/qc-report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID","Raw Reads","Clean Reads","Raw Base","Clean Base","Raw Q30","Clean Q30","Raw GC(%)","Clean GC(%)","Raw Dup(%)","Clean Dup(%)","Adapter(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %stat) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,$stat{$sample}{rawread},$stat{$sample}{cleanread},$stat{$sample}{rawbase},$stat{$sample}{cleanbase},sprintf("%.2f",100*$stat{$sample}{rawq30}),sprintf("%.2f",100*$stat{$sample}{cleanq30}),sprintf("%.2f",100*$stat{$sample}{rawGC}),sprintf("%.2f",100*$stat{$sample}{cleanGC}),$stat{$sample}{rdup},$stat{$sample}{cdup},sprintf("%.2f",100*$stat{$sample}{ada}),$stat{$sample}{pollutions},$stat{$sample}{pollutionr},$stat{$sample}{pollutionn}),"\n";
}
close Out;
open Out,">$dOut/6.report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID","Raw Data(bp)","Clean Data(bp)","Clean Reads(reads)","CleanBasePer(%)","Q30(%)","GC(%)","Dup(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %stat) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,$stat{$sample}{rawbase},$stat{$sample}{cleanbase},$stat{$sample}{cleanread},sprintf("%.2f",100*$stat{$sample}{cleanbase}/$stat{$sample}{rawbase}),sprintf("%.2f",100*$stat{$sample}{cleanq30}),sprintf("%.2f",100*$stat{$sample}{cleanGC}),$stat{$sample}{cdup},$stat{$sample}{pollutions},$stat{$sample}{pollutionr},$stat{$sample}{pollutionn}),"\n";
}
close Out;
my $dir=basename($qc);
mkdir "$dOut/fig" if (!-d "$dOut/fig");
`ln -s $qc/fig/*.png $dOut/fig`;
`ln -s $qc/fig/*.pdf $dOut/fig`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


sub USAGE {
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
		-rawdup		<file>	raw dup list
		-qc	<dir>	fqstp qc dir
		-cleandup	<file>	clean dup list
		-pollution	<file>	pollution summary
		-o	<dir>	QC report
USAGE
	print $usage;
	exit;}
