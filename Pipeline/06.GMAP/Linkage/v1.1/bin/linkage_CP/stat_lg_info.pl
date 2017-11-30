#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ==============================================================
# Get Options
# ==============================================================
my ($in,$o);

GetOptions(
				"help|?" =>\&USAGE,
				"in:s"=>\$in,
				"o:s"=>\$o,
				) or &USAGE;
&USAGE unless ($in);
#===============================================================
# Default optional value 
#===============================================================
#$o||="./map.summary";
my $fOut = defined $o ? "$o/map.stat.xls" : "./map.stat.xls";
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

my (%maxgap,%gappos,%gapmarker,%nloc,%md,%gap_5);
foreach my $fIn (@infiles) {
	my $gappos=0;
	my $curpos=0;
	my $maxgap=0;
	my $gap_5=0;
	my $formerdis=0;
	my $gapmarker="";
	my $formermarker="";
	open IN,$fIn or die $!;
	while (<IN>) {
		chomp;
		next if ($_=~/group/ or $_=~/^$/);
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
	}
	close IN;
	my $name=basename($fIn);
	$name=~s/\.map//g;
	$name=~s/(.*?)\.(\w+)//g;
	$maxgap{$1}{$2}=$maxgap;
	$gappos{$1}{$2}=$gappos;
	$gapmarker{$1}{$2}=$gapmarker;
	$nloc{$1}{$2}=$curpos;
	$md{$1}{$2}=$formerdis;
	$gap_5{$1}{$2}=$gap_5;
}
open OUT,">$fOut" or die $!;
print OUT "#LG\tsex\tnloc\tgap_5\tgap_5(%)\tdis\taver.dis\tmaxgap\tgappos\tgapmarker\n";
my @type=("male","female","sexAver");
foreach my $sex (sort{$a cmp $b} @type) {
	foreach my $lg (sort keys %maxgap) {
		my ($gap_5_ra,$aver_dis);
		if ($nloc{$lg}{$sex}>0) {
			$gap_5_ra = sprintf("%.3f", $gap_5{$lg}{$sex}/$nloc{$lg}{$sex});
			$aver_dis = sprintf("%.3f", $md{$lg}{$sex}/($nloc{$lg}{$sex}-1));
		}else{
			$gap_5_ra =0;$aver_dis = 0;
		}
		print OUT "$lg\t$sex\t$nloc{$lg}{$sex}\t$gap_5{$lg}{$sex}\t$gap_5_ra\t$md{$lg}{$sex}\t$aver_dis\t$maxgap{$lg}{$sex}\t$gappos{$lg}{$sex}\t$gapmarker{$lg}{$sex}\n";
	}
}
print OUT "\tnloc_s\tgap_5\tmd\tnloc_f\tgap_5\tmd\tnloc_m\tgap_5\tmd\tmaxgap\n";
print OUT "Total\t".sum(\%nloc,'sexAver')."\t".sum(\%gap_5,'sexAver')."\t".sum(\%md,'sexAver')."\t".sum(\%nloc,'female')."\t".sum(\%gap_5,'female')."\t".sum(\%md,'female')."\t".sum(\%nloc,'male')."\t".sum(\%gap_5,'male')."\t".sum(\%md,'male')."\t";
print OUT max(\%maxgap,'sexAver')."\n";
close OUT;

# ==============================================================
# sub function
# ==============================================================
sub max{
my ($hash,$ele)=@_;
my $max=0;
foreach my $key(sort keys %$hash){
    $max=$$hash{$key}{$ele} if($$hash{$key}{$ele}>$max);
}
return $max;
}
sub sum{
my ($hash,$ele)=@_;
my $sum=0;
foreach my $key(sort keys %$hash){
    $sum+=$$hash{$key}{$ele} if(exists $$hash{$key}{$ele});
}
return $sum;
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Shi Tongwei <shitw\@biomarker.com.cn> 
=====================================================================================================
Discription:
	This script is used to  calculate max gap in a map
=====================================================================================================

Usage:
  Options:
  -i	<file>	required	input map file or directory
  -o	<str>	optional	output file ,default [./map.stat.xls]
  -h		Help

USAGE
	print $usage;
	exit;
}

