#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gff , $anno , $fOut ,$fregion);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"region:s"=>\$fregion,
	"anno:s"=>\$anno,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fregion and $gff and $anno and $fOut );
open In,$fregion;
my %region;
my $region;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	if (/^@/) {
		my ($chr,$pos1,$pos2)=split(/\s+/,$_);
		$region=join("\t",$pos1,$pos2);
	}else{
		$region{$chr}{$region}{eff}{$_}=1;
	}
}
close In;
open In,$anno;
open Out,">$fOut.summary";
my %kdetail;
my %gdetail;
my %enrich;
my %edetail;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($id,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$eid,$eanno)=split(/\t/,$_);
	my ($id,$chr,$pos1,$pos2)=split(/\:/,$id);
	my $regioned=0;
	foreach my $region (sort keys $region{$chr}) {
		my ($pos3,$pos4)=split(/\t/,$region);
		if (($pos1 > $pos3 && $pos1 < $pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos1 > $pos3 && $pos2 < $pos4)) {
			$regioned=1;
		}
	}
	my $info=$info{$id};
	my @koid=split(/:/,$koid);
	my @kdetail=split(/:/,$koanno);
	for (my $i=0;$i<@koid;$i++) {
		next if ($koid[$i] eq "--");
		$enrich{$koid[$i]}{total}++;
		$enrich{$koid[$i]}{enrich}++ if($regioned ==1);
		$kdetail{$koid[$i]}=$kdetail[$i];
	}
	my @goid=split(/,/,$goid);
	my @gdetail=split(/:/,$goanno);
	for (my $i=0;$i<@goid;$i++) {
		next if ($goid[$i] eq "--");
		$enrich{$goid[$i]}{total}++;
		$enrich{$goid[$i]}{enrich}++ if($regioned ==1);
		$gdetail{$goid[$i]}=$gdetail[$i];
	}
	my @eid=split(/,/,$eid);
	my @edetail=split(/;/,$eanno);
	for (my $i=0;$i<@eid;$i++) {
		next if ($eid[$i] eq "--");
		$enrich{$eid[$i]}{total}++;
		$enrich{$eid[$i]}{enrich}++ if($regioned ==1);
		$edetail{$eid[$i]}=$edetail[$i];
	}
}
close Out;
close In;
open Out,">$fOut.kegg.stat";
foreach my $koid (sort keys %kdetail) {
	$enrich{$koid}{enrich}||=0;
	$enrich{$koid}{total}||=0;
	print Out join("\t",$koid,$kdetail{$koid},$enrich{$koid}{enrich},$enrich{$koid}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;
open Out,">$fOut.go.stat";
foreach my $goid (sort keys %gdetail) {
	$enrich{$goid}{enrich}||=0;
	$enrich{$goid}{total}||=0;
	print Out join("\t",$goid,$gdetail{$goid},$enrich{$goid}{enrich},$enrich{$goid}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;
open Out,">$fOut.eggnog.stat";
foreach my $eggnog (sort keys %edetail) {
	$enrich{$eggnog}{enrich}||=0;
	$enrich{$eggnog}{total}||=0;
	print Out join("\t",$eggnog,$edetail{$eggnog},$enrich{$eggnog}{enrich},$enrich{$eggnog}{total},scalar keys %eff,scalar keys %stat),"\n";
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -snp	<file>	input snp file name
  -indel	<file>	input indel file name
  -anno	<file>	input anno file name
  -out	<key>	output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}
