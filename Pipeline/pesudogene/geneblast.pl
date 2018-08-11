#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($input,$output,$geneblast,$ref);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$ref,
	"geneblast:s"=>\$geneblast,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($ref and $geneblast and $output);
open In,$ref;
$/=">";
my %seq;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@line)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	$seq{$id}=join("",@line);
}
close In;
$/="\n";
open In,$geneblast;
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/([^|]*)\|([^:]*):(\d+)\.\.(\d+)\|([^|])\|([^|]*)\|([^|]*)\|rank:1/) {
		my $geneid=$1;
		my $chrid=$2;
		my $start=$3;
		my $end=$4;
		my $flag=$5;
		my $coverage=$6;
		my $score=$7;
		my $seq=substr($seq{$chrid},$start-1,$end-$start+1);
		if ($flag eq "-") {
			$seq=reverse($seq);
			$seq=~tr/atcgATCG/TAGCTAGC/;
		}
		if (!exists $filehand{$geneid}) {
			open $filehand{$geneid},">$output/$geneid.dna.fa";
		}
		print {$filehand{$geneid}} ">$_\n$seq\n";
	}else{
		next;
	}
}
close In;
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
  -input	<file>	input dir
  -output	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
