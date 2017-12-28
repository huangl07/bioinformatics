#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"g:s"=>\$pop,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $pop);
open In,$pop;
my %group;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$groupid)=split(/\t/,$_);
	$group{$id}=$groupid;
}
close In;
open In,$fIn;
my @indi;
my %out;
my @head;
push @head,"ind";
push @head,"pop";
my $nloc=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	if (/^#/) {
		push @indi,@geno;	
		for (my $i=0;$i<@geno;$i++) {
			push @{$out{$indi[$i]}},$geno[$i];
			push @{$out{$indi[$i]}},$group{$geno[$i]};
		}
	}else{
		$nloc++;
		push @head,join("_",$chr,$pos);
		push @head,"";
		my @ale=split(",",join(",",$ref,$alt));
		my @format=split(/\:/,$format);
		my %geno;
		for (my $i=0;$i<@ale;$i++) {
			$geno{$i}=100+$i*10;
		}
		for (my $i=0;$i<@geno;$i++) {
			my @info=split(/\:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "GT") {
					if ($info[$j] eq "./.") {
						push @{$out{$indi[$i]}},join(",",0,0);
					}else{
						my @gt=split(/\//,$info[$j]);
						push @{$out{$indi[$i]}},join(",",$geno{$gt[0]},$geno{$gt[1]});
					}
				}
			}
		}
	}
}
close In;
open Out,">$fOut";
print Out join(",",$nloc,scalar keys %out,1),"\n";
print Out "," x scalar @head,"\n";
print Out join(",",@head),"\n";
foreach my $ind (sort keys %out) {
	print Out join(",",@{$out{$ind}}),"\n";
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
