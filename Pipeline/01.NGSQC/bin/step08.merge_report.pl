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
my ($rawlist , $rawdup , $cleanlist , $cleanfig , $cleandup , $dOut,$pollution);
GetOptions(
				"help|?" =>\&USAGE,
				"rawlist:s"=>\$rawlist,,
				"rawdup:s"=>\$rawdup,
				"cleanlist:s"=>\$cleanlist,
				"cleanfig:s"=>\$cleanfig,
				"cleandup:s"=>\$cleandup,
				"pollution:s"=>\$pollution,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($rawlist and $rawdup and $cleanlist and $cleanfig and $cleandup and $dOut);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
my %sample;
open In,$rawlist;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$stat)=split(/\s+/,$_);
	open Stat,$stat;
	while (<Stat>) {
		next if ($_ eq "" ||/^$/||!/Total/);
		(undef,$sample{raw}{$sample},$sample{raw}{$sample}{read},$sample{raw}{$sample}{base},$sample{raw}{$sample}{A},$sample{raw}{$sample}{T},$sample{raw}{$sample}{G},$sample{raw}{$sample}{C},$sample{raw}{$sample}{N},$sample{raw}{$sample}{GC},$sample{raw}{$sample}{Q30},$sample{raw}{$sample}{Q20},$sample{raw}{$sample}{AQ},$sample{raw}{$sample}{enzyme1},$sample{raw}{$sample}{enzyme2})=split(/\s+/,$_);
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
			$sample{raw}{$sample}{dup}=$rdup;
		}
	}
	close Dup;
}
close In;
open In,$cleanlist;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$stat)=split(/\s+/,$_);
	open Stat,$stat;
	while (<Stat>) {
		next if ($_ eq "" ||/^$/||!/Total/);
		(undef,undef,$sample{clean}{$sample}{read},$sample{clean}{$sample}{base},$sample{clean}{$sample}{A},$sample{clean}{$sample}{T},$sample{clean}{$sample}{G},$sample{clean}{$sample}{C},$sample{clean}{$sample}{N},$sample{clean}{$sample}{GC},$sample{clean}{$sample}{Q30},$sample{clean}{$sample}{Q20},$sample{clean}{$sample}{AQ},$sample{clean}{$sample}{enzyme1},$sample{clean}{$sample}{enzyme2})=split(/\s+/,$_);
	}
	close Stat;
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
			$sample{clean}{$sample}{dup}=(split(/\s+/,$_))[2];
		}
	}
	close Dup;
}
close In;
open In,"$pollution";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$po1,$po2,$po3)=split(/\s+/,$_);
	$sample{raw}{$sample}{pollutions}=$po1;
	$sample{raw}{$sample}{pollutionr}=$po2;
	$sample{raw}{$sample}{pollutionn}=$po3;
}
close In;
open Out,">$dOut/qc-report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID\t","Raw Reads","Clean Reads","Raw Base","Clean Base","Raw Q30","Clean Q30","Raw AQ(%)","Clean AQ(%)","Raw GC(%)","Clean GC(%)","Raw Dup(%)","Clean Dup(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %{$sample{raw}}) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,$sample{raw}{$sample}{read},$sample{clean}{$sample}{read},$sample{raw}{$sample}{base},$sample{clean}{$sample}{base},$sample{raw}{$sample}{Q30},$sample{clean}{$sample}{Q30},$sample{raw}{$sample}{AQ},$sample{clean}{$sample}{AQ},$sample{raw}{$sample}{GC},$sample{clean}{$sample}{GC},$sample{raw}{$sample}{pollutions},$sample{raw}{$sample}{pollutionr},$sample{raw}{$sample}{pollutionn}),"\n";
}
close Out;
open Out,">$dOut/10.report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID\t","Raw Base","Clean Base","Clean Reads","CleanBasePer(%)","Q30(%)","GC(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %{$sample{raw}}) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,sprintf("%.2f",$sample{raw}{$sample}{base}/1000000000),sprintf("%.2f",$sample{clean}{$sample}{base}/1000000000),sprintf("%.2f",$sample{clean}{$sample}{read}/1000000),sprintf("%.2f",100*$sample{clean}{$sample}{base}/$sample{raw}{$sample}{base}),$sample{clean}{$sample}{Q30},$sample{clean}{$sample}{GC},$sample{raw}{$sample}{pollutions},$sample{raw}{$sample}{pollutionr},$sample{raw}{$sample}{pollutionn}),"\n";
}
close Out;
my $dir=basename($cleanlist);
mkdir "$dOut/fig" if (!-d "$dOut/fig");
`ln -s $cleanfig/*.png $dOut/fig`;
`ln -s $cleanfig/*.pdf $dOut/fig`;

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
		-rawlist	<file>	raw qc list
		-rawdup		<file>	raw dup list
		-cleanlist	<file>	clean qc list
		-cleanfig	<dir>	clean qc fig
		-cleandup	<file>	clean dup list
		-pollution	<file>	pollution summary
		-o	<dir>	QC report
USAGE
	print $usage;
	exit;}
