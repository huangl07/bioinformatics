#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($slist,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$slist,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($slist and $fOut );
open In,$slist;
my %catalog;
my %depth;
my %Total;
my $TotalDep=0;
while (<In>) {
	next if ($_ eq "" ||/^$/);
	my ($sample,$slist)=split(/\s+/,$_);
	open Slist,"gunzip -c $slist|";
	while (<Slist>) {
		chomp;
		next if ($_ eq "" ||/^$/ || /^#/);
		my (undef,undef,$cata,undef,undef,undef,$dep,undef)=split;
		$catalog{$sample}{$cata}=1;
		$depth{$sample}+=$dep;
		$Total{$cata}=1;
		$TotalDep+=$dep;
	}
	close Slist;
}
close In;
open Out,">$fOut";
print Out "#sample\ttags Number\tAverage Number\n";
foreach my $sample (sort keys %catalog) {
	my @out;
	push @out,$sample;
	push @out,scalar keys %{$catalog{$sample}};
	push @out,sprintf("%.4f",$depth{$sample}/scalar keys %{$catalog{$sample}});
	print Out join("\t",@out),"\n";
}
print Out join("\t","Total",scalar keys %Total,$TotalDep/scalar keys %Total),"\n";
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
  -o	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}
