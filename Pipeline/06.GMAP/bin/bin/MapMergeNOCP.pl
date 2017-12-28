#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$dOut,$adjust);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"out:s"=>\$dOut,
	"adjust"=>\$adjust,
			) or &USAGE;
&USAGE unless ($dmap and $dOut );
mkdir $dOut if (!-d $dOut);
my @map=glob("$dmap/*.out");
open Out,">$dOut/total.map";
my %Marker;
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	my @order;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\t/,$_);
		$Marker{$id}=join(",",$lgID,$pos);
		push @order,$id;
		if ($max < $pos) {
			$max=$pos;
		}
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(80)+120;
	}
	foreach my $id (@order) {
		my ($lgid,$pos)=split(/\,/,$Marker{$id});
		$pos=$newdis/$max*$pos;
		$Marker{$id}=join(",",$lgid,$pos);
		print Out $id,"\t",$pos,"\n";
	}	
}
close Out;
my @marker=glob("$dmap/*.marker");
open Out,">$dOut/total.marker";
open CSV,">$dOut/total.csv";
my $head;
my $chead;
my @out;
my @cout;
foreach my $marker (@marker) {
	open In,$marker;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my @info=split;
		next if (scalar @info < 3);
		if (/MarkerID/) {
			$head=$_;
			my (undef,$nhead)=split(/\s+/,$_,2);
			$nhead=~s/\t/,/g;
			my @head=split(/\,/,$nhead);
			$chead="Genotype,,,".join(",",@head);
		}else{
			push @out,$_;
			my ($id,$info)=split(/\s+/,$_,2);
			$info=~s/\t/,/g;
			$info=~s/X/H/g;
			$info=~s/U/-/g;
			push @cout,join(",",$id,$Marker{$id},$info);

		}
	}
	close In;
}
print Out join("\n",$head,@out);
close Out;
print CSV join("\n",$chead,@cout);
close CSV;

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
  -dmap	<file>	input file name
  -out	<file>	output file
  -adjust	<file>	adjust map distcance
  -h         Help

USAGE
        print $usage;
        exit;
}
