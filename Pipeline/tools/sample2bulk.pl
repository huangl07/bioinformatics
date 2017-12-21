#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$groid);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"g:s"=>\$groid,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $groid);
open In,$groid;
my %groid;
my %id;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$gid)=split(/\s+/,$_);
	$groid{$id}=$gid;
	$id{$gid}=1;
}
close In;
open In,$fIn;
open Out,">$fOut";
my @Indi;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ );
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	if (/^##/) {
		print Out "$_\n";
	}elsif (/^#/) {
		push @Indi,@geno;
		print Out join("\t",$chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,"GT:AD:DP",sort keys %id),"\n";
		next;
	}else{
		my %newgeno;
		my @format=split(/:/,$format);
		for (my $i=0;$i<@geno;$i++) {
			my $nid=$groid{$Indi[$i]};
			my ($gt,$ad,$dp);
			my @info=split(/\:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "GT") {
					$gt=$info[$j];
					next if($gt eq "./.");
					my @gt=split(/\//,$gt);
					$newgeno{$nid}{$gt[0]}++;
					$newgeno{$nid}{$gt[1]}++;
				}
			}
		}
		my @out;
		push @out,join("\t",$chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,"GT:AD:DP");
		my $n=scalar split(/\,/,join(",",$ref,$alt));
		foreach my $id (sort keys %id) {
			my ($gt,$ad,$dp);
			my @gt=sort keys %{$newgeno{$id}};
			if (scalar @gt == 1) {
				$gt="$gt[0]/$gt[0]";
				my @ad;
				for (my $i=0;$i<$n;$i++) {
					if ($i eq $gt[0]) {
						push @ad,$newgeno{$id}{$gt[0]};
					}else{
						push @ad,0;
					}
				}
				$ad=join(",",@ad);
				$dp=$newgeno{$id}{$gt[0]};
			}elsif (scalar @gt==2) {
				$gt="$gt[0]/$gt[1]";
				my @ad;
				for (my $i=0;$i<$n;$i++) {
					if ($i eq $gt[0]) {
						push @ad,$newgeno{$id}{$gt[0]};
					}elsif ($i eq $gt[1]) {
						push @ad,$newgeno{$id}{$gt[1]};
					}else{
						push @ad,0;
					}
				}
				$ad=join(",",@ad);
				$dp=$newgeno{$id}{$gt[0]}+$newgeno{$id}{$gt[1]};
			}else{
				$gt="./.";
				$ad=".";
				$dp=".";
			}
			push @out,join(":",$gt,$ad,$dp);
		}
		print Out join("\t",@out),"\n";
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
