#!/usr/bin/perl -w
# Writer:          Shi tongwei <shitw@biomarker.com.cn>
# Modifier:        qix <qix@biomarker.com.cn>
# Last Modified:   2015/10/15
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.1";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($in,$o,$loc,$cycle);

GetOptions(
				"help|?" =>\&USAGE,
				"c:s"=>\$cycle,
				"in:s"=>\$in,
				"loc:s"=>\$loc,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($in);
#===============================================================
# Default optional value 
#===============================================================
my $fOut = defined $o ? "$o/MAP.STAT.xls" : "./MAP.STAT.xls";
my @infiles;
if (-d $in) {
	@infiles=glob("$in/*.map");
	if (scalar @infiles == 0) {
		@infiles=glob("$in/*.map");
	}
}elsif(-f $in){
	$infiles[0]=$in;
}else{
die "no files input!!!\n" ;
}
open LOC,$loc or die $!;
my ($ind,%miss);
while (<LOC>) {
	chomp;
	next if (/^\s*$/||$_ eq ""||/=/||/;/);
	my ($marker,$type,$type2,@genotype) = split /\s+/ , $_ ;
	$ind = scalar(@genotype);
	my $nmiss=0;
	for ( my $i=0;$i<@genotype ;$i++) {
		if ($genotype[$i] eq "--") {
			$nmiss++;
		} 
	}
	$miss{$marker}+=$nmiss;
}
#print Dumper %miss;die;
close LOC;
my (%maxgap,%gappos,%gapmarker,%gap_5,%nloc,%md,%summiss);  #modify by qix ,add %gap_5
foreach my $fIn (@infiles) {
	my $gappos=0;
	my $curpos=0;
	my $maxgap=0;
	my $gap_5=0;                                   
	my $formerdis=0;
	my $gapmarker="";
	my $formermarker="";
	my $sexmiss=0;
	my @sexmarker=();
	open IN,$fIn or die $!;
	while (<IN>) {
		chomp;
		next if ($_=~/group/ or $_=~/^$/); #已经包含连锁群为空的情况
		my ($marker,$dis)=split /\s+/,$_;
		my $gap=$dis-$formerdis;
		if ($gap < 5) {                            
			$gap_5 ++ ;
		}
		$curpos++;
		if ($gap>$maxgap) {
			$maxgap=$gap;
			$gappos=$curpos;
			$gapmarker=$formermarker."_".$marker;
		}
		$formermarker=$marker;
		$formerdis=$dis;
		push @sexmarker ,$marker;
	}	
	close IN;
	foreach my $key (@sexmarker) { #找缺失MARKER
		if (exists $miss{$key}) {
			$sexmiss+=$miss{$key};
		}
	}
#	print $sexmiss;die;
	my $name=basename($fIn);
	$name=~s/\.map//g;
	$name=~s/(.*?)\.(\w+)//g;
#	print $1,$2;die; $1为连锁群号，$2为性别
	$maxgap{$2}=$maxgap;
	$gappos{$2}=$gappos;
	$gapmarker{$2}=$gapmarker;
	$nloc{$2}=$curpos;
	$md{$2}=$formerdis;
	$gap_5{$2}=$gap_5;                     #modify by qix ,add %gap_5
	$summiss{$2}=$sexmiss;
}
#print Dumper %maxgap;die;

open OUT,">>$fOut" or die $!;
foreach my $sex (sort{$a cmp $b}keys %maxgap) {	
	my ($sumnum,$nmiss_ra,$gap_5_ra,$aver_dis); 
	if ($nloc{$sex}>0) {
		$sumnum =  $nloc{$sex}*$ind;
		$nmiss_ra = sprintf "%.3f", $summiss{$sex}/($nloc{$sex}*$ind);
		$gap_5_ra = sprintf "%.3f", $gap_5{$sex}/$nloc{$sex};
		$aver_dis = sprintf "%.3f", $md{$sex}/($nloc{$sex}-1);
	}else{
		$sumnum = 0;$nmiss_ra =0;$gap_5_ra =0;$aver_dis = 0;
	}
	print OUT "cycle$cycle\t$sex\t$nloc{$sex}\t$ind\t$sumnum\t$summiss{$sex}\t$nmiss_ra\t$gap_5{$sex}\t$gap_5_ra\t$md{$sex}\t$aver_dis\t$maxgap{$sex}\t$gappos{$sex}\t$gapmarker{$sex}\n";
}
close OUT;

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact:	Qixin <qix\@biomarker.com.cn> 
=====================================================================================================
Discription:
	This script is used to  calculate max gap in a map
=====================================================================================================

Usage:
  Options:
  -i	<file>	required	input map file or directory
  -o	<str>	optional	output file ,default [./MAP.STAT.xls]
  -h		Help

USAGE
	print $usage;
	exit;
}

