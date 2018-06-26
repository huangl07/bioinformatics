#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($in,$out,$list);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$in,
	"o:s"=>\$out,
	"list:s"=>\$list,
			) or &USAGE;
&USAGE unless ($in);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        select read from bam file
        eg:
        perl $Script -i  -o -list

Usage:
  Options:
  -i	<file>	input fq file
  -o	<file>	output result file
  -list	<str>	fq's id list 
  -h		Help

USAGE
        print $usage;
        exit;
}

my %key;
open IN,"$in";
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($ids,@others)=split/\t/,$_;
    $key{$ids}=1;
    }
close IN;

open In,$in;
if($in=~/gz$/){
	open In,"zcat  $in|";
}
$/="\@";
open Out,">$out";
while(<In>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/);
	my $seq=$_;
	my $seqid=(split(/\s+/,(split(/\n/,$seq))[0]))[0];
    if(exists $key{$seqid}){
        print Out "\@$seq";
        }
}
close In;
close Out;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
