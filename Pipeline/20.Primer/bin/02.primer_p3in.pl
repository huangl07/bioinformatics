#!/mnt/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "input:s"=>\$ref,
        "output:s"=>\$out,
        "type:s"=>\$type,
                        ) or &USAGE;
&USAGE unless ($ref and $out and $type);
my %stat;
my $number;
open Misa,"$out.$type.misa";
#open Out,">>$out/$type.new.fa";
while (<Misa>){
	chomp;
	next if ($_ eq "" || /^$/|| /^#/);#chr:ID	snp.nr	type	snp	size	start	end
	my ($ids,$nr,$number,$marker,$ref,$alt,$size,$start,$end,$nstart,$nend)=split(/\t/,$_);
	$stat{$number}{chr}=(split(/\:/,$ids))[0];
	$stat{$number}{chrID}=join("",$ids,"_",$nr,"_",$ref);
	$stat{$number}{alt}=$alt;
	$stat{$number}{size}=$size;
	$stat{$number}{nstart}=$nstart;
	$stat{$number}{nend}=$nend;
}
close Misa;
open Out,">$out.$type.new.fa";
open P3in,">$out.$type.p3in";
open In,$ref;
my ($id,$seq);
while (<In>){
	chomp;
	next if ($_ eq "" || /^$/|| /^#/);
	$id=$_ ;
	$seq=<In>;
	$id =(split(/\>/,$id))[1];
	my $len = length$seq;
	foreach $number (sort { $a <=> $b } keys %stat){
		if ($stat{$number}{chr} eq $id){
			my $rend;
			if ($type=~/snp/){
				$rend=1001;
			}else{
				$rend =$stat{$number}{size} + 300 ;
			}
			my $result  = substr($seq,$stat{$number}{nstart},$rend);
			print Out $stat{$number}{chrID},"\n";
			print Out $result,"\n";
			if ($type=~/snp/){
				print P3in "PRIMER_SEQUENCE_ID=$stat{$number}{chrID}\nSEQUENCE_TEMPLATE=$result\n";
				print P3in "PRIMER_PRODUCT_SIZE_RANGE=600-800\n";
				print P3in "TARGET=",450,",",100,"\n";
				print P3in "PRIMER_MAX_END_STABILITY=780\n=\n";
			}else{
				print P3in "PRIMER_SEQUENCE_ID=$stat{$number}{chrID}\nSEQUENCE_TEMPLATE=$result\n";
				print P3in "PRIMER_PRODUCT_SIZE_RANGE=100-200\n";
				print P3in "TARGET=",80,",",$stat{$number}{size}+70,"\n";
				print P3in "PRIMER_MAX_END_STABILITY=150\n=\n";}
		}
	}
}
close In;
close Out;
close P3in;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
########################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c

Usage:
  Options:
	-input	<file>  input file name
	-output	<file>  split windows sh
	-type	<str>	marker type
	-h			Help

USAGE
        print $usage;
        exit;
}
	

