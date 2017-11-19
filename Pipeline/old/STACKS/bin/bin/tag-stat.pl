#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ulist,$group,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"u:s"=>\$ulist,
	"g:s"=>\$grouplist,
	"o:s"=>\$sample,
			) or &USAGE;
&USAGE unless ($ulist and $group);
my %group;
if ($grouplist) {
	open In,$grouplist;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$group)=split(/\s+/,$_);
		my @group=split(/\,/,$group);
		foreach my $sample (@group) {
			$group{$sample}=$group;
		}
	}
	close In;
}
open In,$ulist;
my %dep;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$ulist)=split(/\s+/,$_);
	my $gID;
	if ($grouplist) {
		$gID=$group{$sample};
	}else{
		$gID=1;
	}
	open TAG,"gunzip -c $ulist.tags.tsv.gz|";
	my $total=0;
	my $dep=0;
	while (<TAG>) {
		chomp;
		next if ($_ eq ""||/^$/);
		if (/consensus/) {
			$total++;
		}
		if (/primary/ || /secondary/) {
			$ndp++;
		}
	}
	close In;
	$dep{$gID}{$sample}=$total/$ndp;
}
close In;
open Out,">$fOut";
foreach my $group (sort keys %dep) {
	my @sample=sort {$dep{$group}{$a}<=$dep{$group}{$b}} keys %{$dep{$group}};
	if (scalar @sample >1) {
		print Out join("\n",$sample[0],$sample[-1]),"\n";
	}else{
		print Out join("\n",$sample[0])
	}
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	group list output file
  -g	<file>	group list
  -h         Help

USAGE
        print $usage;
        exit;
}
