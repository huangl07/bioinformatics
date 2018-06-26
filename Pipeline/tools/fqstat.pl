#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output
			) or &USAGE;
&USAGE unless ($input and $output);
my @fq=glob("$input/*.fastq");
my %region;
foreach my $fq (@fq) {
	my $fqname=basename($fq);
	$fqname=~s/\.fastq//;
	open In,$fq;
	my $seq;
	my $qual;
	while ($seq=<In>) {
		$seq=<In>;
		$qual=<In>;
		$qual=<In>;
		chomp $seq;
		chomp $qual;
		$region{$fqname}{readnum}++;
		$region{$fqname}{basenum}+=length($seq);
		my @qual=split(//,$qual);
		for (my $i=0;$i<=@qual;$i++) {
			my $qual = ord($qual[$i])-33;
			$region{$fqname}{q30}++ if($qual >30);
		}
	}
	close In;
}
open Out,">$output";
print Out "#region\treadnum\tbasenum\tq30\n";
foreach my $region (sort keys %region) {
	print Out join("\t",$region,$region{$region}{readnum},$region{$region}{basenum},$region{$region}{q30}/$region{$region}{basenum}),"\n";
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
  -input	<dir>	input A file name
  -output	<file>	input B file name
  -h         Help

USAGE
        print $usage;
        exit;
}
