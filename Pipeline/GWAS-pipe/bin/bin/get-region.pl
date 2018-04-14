#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut and $block);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"b:s"=>\$block,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $block );
open In,$block;
my %region;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my (undef,@marker)=split(/\s+/,$_);
	my @pos;
	my $chrom ="";
	for (my $i=0;$i<@marker;$i++) {
		my ($chr,$pos)=split(/\_/,$marker[$i]);
		push @pos,$pos;
		die "error block" if ($chrom ne $chr && $chrom ne "");
		$chrom=$chr;
	}
	@pos=sort{$a<=>$b}@pos;
	for (my $i=0;$i<@marker;$i++) {
		$region{$marker[$i]}=join("\t",$chr,$pos[0],$pos[1]);
	}
}
close In;
open In,$fIn;
open Out,">$fOut";
my %snp;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($taxa,$chr,$pos,$maf,$pvalue)=split(/\s+/,$_);
	my $id=join("_",$chr,$pos);
	if (exists $region{$marker[$i]}) {
		print Out $region{$marker[$i]},"\t","blocked\n";
	}else{
		print Out $chr,$pos-500000,$pos+500000,"\t","unblocked\n";
	}
}
close In;
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
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
