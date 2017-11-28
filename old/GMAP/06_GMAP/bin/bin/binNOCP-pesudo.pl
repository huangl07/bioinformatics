#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fgeno,$dir,$Key,$fpos,$winsize,$stepsize);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fgeno,
	"p:s"=>\$fpos,
	"o:s"=>\$dir,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($fgeno and $dir and $Key);
mkdir $dir if (!-d $dir);
open In,$fpos;
my %pos;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/||/^$/);
	my ($id,$chr,$pos)=split(/\s+/,$_);
	$pos{$id}=join("\t",$id,$chr,$pos);
}
close In;
open In,$fgeno;
open Male,">$dir/$Key.male.matrix";
open Mpos,">$dir/$Key.male.pos";
open Female,">$dir/$Key.female.matrix";
open Fpos,">$dir/$Key.female.pos";
my @Head;
while (<In>) {
	#进行每个marker的来源判断，0代表无法确定，1代表第二条，-1代表第1条染色体
	#同时进行拟测交处理，未判断连锁相，maybe not right
	chomp;
	next if ($_ eq "" || /^$/  );
	if (/^#/) {
		(undef,undef,@Head)=split(/\s+/,$_);
		print Male $_,"\n";
		print Female $_,"\n";
		print Mpos "#MarkerID\tChr\tPos\n";
		print Fpos "#MarkerID\tChr\tPos\n";
	}else{
		my ($id,$type,@geno)=split(/\s+/,$_);
		my @male;
		my @female;
		if ($type eq "aaxbb") {
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "aa") {#nn或ll
					push @male,"aa 1";
					push @female,"aa 1";
				}elsif ($geno[$i] eq "bb") {#np或ll
					push @male,"bb -1";
					push @female,"bb -1";
				}elsif ($geno[$i] eq "ab") {#lm或nn
					push @male,"ab 1";
					push @female,"ab -1";
				}else{
					push @male,"-- 0";
					push @female,"-- 0";
				}
			}
			print Male "$id\t$type\t",join("\t",@male),"\n";
			print Mpos $pos{$id},"\n";
			print Female "$id\t$type\t",join("\t",@female),"\n";
			print Fpos $pos{$id},"\n";
		}
	}
}
close In;
close Male;
close Female;
close Mpos;
close Fpos;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description: pesudo cross
	eg:
	perl $Script -i -o -k

Usage:
  Options:
  Options:
  -i	<file>	input genotype file name
			#MakrerID\tTYPE\tSample
  -p	<file>	input pos file
			#MakrerID\tChr\tPos
  -o	<dir>	output dir
  -k	<str>	output file name 
  -h         Help
  
  -h         Help

USAGE
        print $usage;
        exit;
}