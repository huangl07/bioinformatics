#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fanno,$statfile,$Key);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fanno,
	"o:s"=>\$statfile,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($fanno and $statfile and $Key);
open In,$fanno;
my %type;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	my ($chr,$pos,$ref,$alt,$postion,$geneID,$Mut,undef)=split(/\t/,$_);
	my $info=join("\t",$chr,$pos);
	$type{pos}{$postion}{$info}++;
	$type{exonic}{$Mut}{$info}++ if ($Mut ne "-");
}
close In;
close Out;
open Out,">$statfile";
print Out "#pos\t$Key\n";
my $sum=0;
foreach my $type (sort keys %{$type{pos}}) {
	print Out $type,"\t",scalar keys %{$type{pos}{$type}},"\n";
	$sum+=scalar keys %{$type{pos}{$type}};
}
print Out "Total\t$sum\n";
print Out "#mutype\t$Key\n";
$sum=0;
foreach my $type (sort keys %{$type{exonic}}) {
	print Out $type,"\t",scalar keys %{$type{exonic}{$type}},"\n";
	$sum+=scalar keys %{$type{exonic}{$type}};
}
print Out "Total\t$sum\n";
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -i	<file>	input file name
  -o	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
