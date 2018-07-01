#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"p:s"=>\$pop,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
mkdir $fOut if (!-d $fOut);
my %pop;
open In,$pop;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^#/);
	my ($id,$popid)=split(/\s+/,$_);
	$pop{$id}=$popid;
}
close In;
open In,$fIn;
my %seq;
my @Indi;
my %BASE=(
	"AA"=>"A","GG"=>"G","CC"=>"C","TT"=>"T",
	"AT"=>"W","AG"=>"R","AC"=>"M",
	"CG"=>"S","CT"=>"Y",
	"GT"=>"K"
);
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);
		foreach my $indi (@Indi) {
			open $filehand{$indi},">$fOut/$indi.msmc";
		}
	}else{
		my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@indi)=split(/\t/,$_);
		my @ale=split(/\,/,join(",",$REF,$ALT));
		for (my $i=0;$i<@indi;$i++) {
			my $geno=(split(/\:/,$indi[$i]))[0];
			my ($g1,$g2)=split(/\//,$geno);
			if ($g1 eq ".") {
				print {$filehand{$Indi[$i]}} join("\t",$CHROM,$POS,1,"?"),"\n";
			}else{
				if (!exists $BASE{join("",sort($ale[$g1],$ale[$g2]))}){
					print {$filehand{$Indi[$i]}} join("\t",$CHROM,$POS,1,"?"),"\n";
				}else{
					print {$filehand{$Indi[$i]}} join("\t",$CHROM,$POS,1,$BASE{join("",sort($ale[$g1],$ale[$g2]))}),"\n";
				}
			}
		}
	}
}
close In;
foreach my $indi (@Indi) {
	close $filehand{$indi};
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
