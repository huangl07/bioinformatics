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
open Other,">$dir/$Key.other.matrix";
open Opos,">$dir/$Key.other.pos";
my @Head;
while (<In>) {
	#进行每个marker的来源判断，0代表无法确定，1代表第二条，-1代表第1条染色体
	#同时进行拟测交处理，未判断连锁相，maybe not right
	chomp;
	next if ($_ eq "" || /^$/ || /aaxbb/ );
	if (/^#/) {
		(undef,undef,@Head)=split(/\s+/,$_);
		print Other $_,"\n";
		print Male $_,"\n";
		print Female $_,"\n";
		print Mpos "#MarkerID\tChr\tPos\n";
		print Fpos "#MarkerID\tChr\tPos\n";
		print Opos "#MarkerID\tChr\tPos\n";
	}else{
		my ($id,$type,@geno)=split(/\s+/,$_);
		my @male;
		my @female;
		my @other;
		my $miss=0;
		if ($type eq "abxcd") {#abxcd,b为m;d为p
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "ac") {#nn或ll
					push @male,"ac -1";
					push @female,"ac -1";
				}elsif ($geno[$i] eq "ad") {#np或ll
					push @male,"ad -1";
					push @female,"ad 1";
				}elsif ($geno[$i] eq "bc") {#lm或nn
					push @male,"bc 1";
					push @female,"bc -1";
				}elsif ($geno[$i] eq "bd") {#np或lm
					push @male,"bc 1";
					push @female,"bc 1";
				}else{
					$miss++;
					push @male,"-- 0";
					push @female,"-- 0";
				}
			}
		}elsif ($type eq "efxeg") {#efxeg,f为m，g为p
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "ee") {#nn,lm
					push @male,"ee -1";
					push @female,"ee -1";
				}elsif ($geno[$i] eq "ef") {#lm
					push @male,"ef 1";
					push @female,"ef -1";
				}elsif ($geno[$i] eq "eg") {#np
					push @male,"bc -1";
					push @female,"bc 1";
				}elsif ($geno[$i] eq "bd") {#lm/np
					push @male,"bc 1";
					push @female,"bc 1";
				}else{
					$miss++;
					push @male,"-- 0";
					push @female,"-- 0";
				}
			}
		}elsif ($type eq "hkxhk") {#hkxhk 不做拟测交处理
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "hh") {
					push @other,"hh 0 0";
				}elsif ($geno[$i] eq "kk") {
					push @other,"kk 0 0";
				}elsif ($geno[$i] eq "hk") {
					push @other,"hk 0 0";
				}else{
					$miss++;
					push @other,"-- 0 0";
				}
			}
		}elsif ($type eq "lmxll") {
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "lm") {
					push @male,"lm 1";
				}elsif ($geno[$i] eq "ll") {
					push @male,"ll -1";
				}else{
					$miss++;
					push @male,"-- 0";
				}
			}
		}elsif ($type eq "nnxnp") {
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "nn") {
					push @female,"nn -1";
				}elsif ($geno[$i] eq "np") {
					push @female,"np 1";
				}else{
					$miss++;
					push @female,"-- 0";
				}
			}
		}elsif ($type eq "abxcc") {#abxcc b为m
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "ac") {
					push @male,"ac -1";
				}elsif ($geno[$i] eq "bc") {
					push @male,"bc 1";
				}else{
					$miss++;
					push @male,"-- 0";
				}
			}
		}elsif ($type eq "ccxab") {
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "ac") {
					push @female,"ac -1";
				}elsif ($geno[$i] eq "bc") {
					push @female,"bc 1";
				}else{
					$miss++;
					push @female,"-- 0";
				}
			}
		}else{
			die "error type: $id\t$type\n";
		}
		if (scalar @male > 0) {
			print Male "$id\t$type\t",join("\t",@male),"\n";
			print Mpos $pos{$id},"\n";

		}
		if (scalar @female > 0) {
			print Female "$id\t$type\t",join("\t",@female),"\n";
			print Fpos $pos{$id},"\n";

		}
		if (scalar @other > 0) {
			print Other "$id\t$type\t",join("\t",@other),"\n";
			print Opos $pos{$id},"\n";
		}
	}
}
close In;
close Male;
close Female;
close Other;
close Mpos;
close Fpos;
close Opos;



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
