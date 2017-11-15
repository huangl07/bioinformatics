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
my ($rawlist , $rawdup , $cleanlist , $cleandup , $dOut,$config,$adaplog,$pollution,$Key);
GetOptions(
				"help|?" =>\&USAGE,
				"rawlist:s"=>\$rawlist,,
				"rawdup:s"=>\$rawdup,
				"cleanlist:s"=>\$cleanlist,
				"cleandup:s"=>\$cleandup,
				"adaplog:s"=>\$adaplog,
				"pollution:s"=>\$pollution,
				"key:s"=>\$Key,
				"out:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($rawlist and $rawdup and $cleanlist and $cleandup and $dOut);
mkdir $dOut if (!-d $dOut);
open In,$rawlist;
my %sample;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$stat)=split(/\s+/,$_);
	open Stat,$stat;
	while (<Stat>) {
		next if ($_ eq "" ||/^$/ ||!/Total/);
		(undef,undef,$sample{raw}{$sample}{read},$sample{raw}{$sample}{base},$sample{raw}{$sample}{A},$sample{raw}{$sample}{T},$sample{raw}{$sample}{G},$sample{raw}{$sample}{C},$sample{raw}{$sample}{N},$sample{raw}{$sample}{GC},$sample{raw}{$sample}{Q30},$sample{raw}{$sample}{Q20},$sample{raw}{$sample}{AQ},$sample{raw}{$sample}{enzyme1},$sample{raw}{$sample}{enzyme2})=split(/\s+/,$_);

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
			$sample{raw}{$sample}{dup}=(split(/\s+/,$_))[2];
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
		next if ($_ eq "" ||/^$/ ||!/Total/);
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
open In,$adaplog;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$rmAda)=split(/\s+/,$_);
	open Ada,$rmAda;
	my $sum=0;
	my $ada=0;
	while (<Ada>) {
		chomp;
		if (/Pairs Processed:/) {
			$sum=(split(/\s+/,$_))[-1];
		}elsif (/Pairs With Adapters:/) {
			$ada=(split(/\s+/,$_))[-1];
		}
	}
	close Ada;
	$sample{raw}{$sample}{Adap}=$ada/$sample{raw}{$sample}{read};#未知原因为0
}
close In;
open In,"$pollution";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$unif)=split(/\s+/,$_,2);
	$sample{raw}{$sample}{pollution}=$unif;
}
close In;
open Out,">$dOut/$Key.qc-report.xls";
print Out join("\t","#LaneID\tLibraryID\tProjectID\tMajorID\tSampleID\t","Raw Reads","Clean Reads","Raw Base","Clean Base","Raw Q30","Clean Q30","Raw AQ(%)","Clean AQ(%)","Raw GC(%)","Clean GC(%)","Raw Dup(%)","Clean Dup(%)","Adapter(%)","Pollution-species","Pollution-rate(%)"),"\n";
foreach my $sample (sort keys %{$sample{raw}}) {
	my $samples=$sample;
	$samples=~s/\:/\t/g;
	print Out join("\t",$samples,$sample{raw}{$sample}{read},$sample{clean}{$sample}{read},$sample{raw}{$sample}{base},$sample{clean}{$sample}{base},$sample{raw}{$sample}{Q30},$sample{clean}{$sample}{Q30},$sample{raw}{$sample}{AQ},$sample{clean}{$sample}{AQ},$sample{raw}{$sample}{GC},$sample{clean}{$sample}{GC},$sample{raw}{$sample}{dup},$sample{clean}{$sample}{dup},sprintf("%.2f",int($sample{raw}{$sample}{Adap}*10000)/100),$sample{raw}{$sample}{pollution}),"\n";
}
close Out;
my $dir=basename($cleanlist);
mkdir "$dOut/fig" if (!-d "$dOut/fig");
`ln -s $dir/*.png $dOut/fig`;
`ln -s $dir/*.pdf $dOut/fig`;

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

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
		-rawlist	<file>	raw qc list
		-rawdup		<file>	raw dup list
		-cleanlist	<file>	clean qc list
		-cleandup	<file>	clean dup list
		-adaplog	<file>	rm ada list
		-pollution	<file>	pollution summary
		-out	<dir>	QC report
USAGE
	print $usage;
	exit;
}
