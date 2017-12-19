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
my %sample;
open In,"$qc/stat.list";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$stat)=split(/\s+/,$_);
	open Stat,$stat;
	while (<Stat>) {
		next if ($_ eq "" ||/^$/||/#/);
		($sample{raw}{$sample},$sample{raw}{$sample}{rawread},$sample{raw}{$sample}{rawbase},$sample{raw}{$sample}{rawq20},$sample{raw}{$sample}{rawq30},$sample{raw}{$sample}{rawGC},$sample{raw}{$sample}{ada},$sample{raw}{$sample}{cleanread},$sample{raw}{$sample}{cleanbase},$sample{raw}{$sample}{cleanq20},$sample{raw}{$sample}{cleanq30},$sample{raw}{$sample}{cleanGC})=split(/\s+/,$_);
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
open In,$pollution;
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
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID\t","Raw Reads","Clean Reads","Raw Base","Clean Base","Raw Q30","Clean Q30","Raw GC(%)","Clean GC(%)","Raw Dup(%)","Clean Dup(%)","Adapter(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %{$sample{raw}}) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,$sample{raw}{$sample}{rawread},$sample{raw}{$sample}{cleanread},$sample{raw}{$sample}{rawbase},$sample{raw}{$sample}{cleanbase},$sample{raw}{$sample}{rawq30},$sample{raw}{$sample}{cleanq30},$sample{raw}{$sample}{rawGC},$sample{raw}{$sample}{cleanGC},$sample{raw}{$sample}{dup},$sample{clean}{$sample}{dup},$sample{raw}{$sample}{ada},$sample{raw}{$sample}{pollutions},$sample{raw}{$sample}{pollutionr},$sample{raw}{$sample}{pollutionn}),"\n";
}
close Out;
open Out,">$dOut/10.report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID\t","Raw Data(G)","Clean Data(G)","Clean Reads(Mreads)","CleanBasePer(%)","Q30(%)","GC(%)","Dup(%)","Pollution-species","Pollution-rate(%)","Pollution-number"),"\n";
foreach my $sample (sort keys %{$sample{raw}}) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,sprintf("%.3f",$sample{raw}{$sample}{rawbase}/1000000000),sprintf("%.3f",$sample{raw}{$sample}{cleanbase}/1000000000),sprintf("%.3f",$sample{raw}{$sample}{cleanread}/1000000),sprintf("%.2f",100*$sample{raw}{$sample}{cleanbase}/$sample{raw}{$sample}{cleanbase}),$sample{raw}{$sample}{cleanq30},$sample{raw}{$sample}{cleanGC},$sample{raw}{$sample}{pollutions},$sample{raw}{$sample}{pollutionr},$sample{raw}{$sample}{pollutionn}),"\n";
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
