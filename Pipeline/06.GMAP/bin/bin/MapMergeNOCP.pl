#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"out:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($dmap and $dOut );
my @map=glob("$dmap/*.out");
open Out,">$dOut/total.map";
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	print Out "group\t$lgID\n";
	open In,$map;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		print Out $_,"\n";
		my ($id,$pos)=split(/\t/,$_);
		$Marker{$id}=join(",",$$lgID,$pos);
	}
	close In;
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
			(undef,$head)=split(/\s+/,$_,2);
			$head=~s/\t/,/g;
			my @head=split(/\,/,$head);
			$chead="Genotype,,,".join(",",@head);
		}else{
			push @out,$_;
			my ($id,$info)=split(/\s+/,$_,2);
			$info=~s/\t/,/g;
			$info=~s/X/H/g;
			$info=~s/U/-/g;
			push @cout join(",",$id,$Marker{$id},$info),"\n"

		}
	}
	close In;
}
print Out join("\n",$head,@out);
close Out;
print Out join("\n",$chead,@cout);
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
  -h         Help

USAGE
        print $usage;
        exit;
}
