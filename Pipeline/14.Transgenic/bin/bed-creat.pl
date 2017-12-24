#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($sv,$out,$bam,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$sv,
	"o:s"=>\$out,
	"r:s"=>\$ref,
			) or &USAGE;
&USAGE unless ($sv and $out and $ref);
$sv=ABSOLUTE_DIR($sv);
$ref=ABSOLUTE_DIR($ref);
open In,"$ref.fai";
my %len;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$len,undef)=split(/\s+/,$_);
	$len{$id}=$len;
}
close In;
open In,$sv;
open Bed,">$out";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	if (/insert/ &&/CTX/) {
		my ($chr1,$pos1,$chr2,$pos2,undef)=split(/\t/,$_);
		my $start;
		my $end;
		if ($chr1 =~ /insert/) {
			$start=$pos2-1000;
			$end=$pos2+1000;
			$start=0 if($start < 0);
			$end=$len{$chr2} if($end > $len{$chr2});
			print Bed join("\t",$chr2,$start,$end),"\n";
			print Bed join("\t",$chr1,0,$len{$chr1});
		}else{
			$start=$pos1-1000;
			$end=$pos1+1000;
			$start=0 if($start < 0);
			$end=$len{$chr1} if($end > $len{$chr1});
			print Bed join("\t",$chr1,$start,$end),"\n";
			print Bed join("\t",$chr2,0,$len{$chr2});
		}
	}
}
close In;
close Bed;

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

  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -pop	<str>	group list
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
