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
my @map=glob("$dmap/*.sexAver.map");
open Out,">$dOut/total.sexAver.map";
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	print Out "group\t$lgID\n";
	open In,$map;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		print Out $_,"\n";
	}
	close In;
}
close Out;
@map=glob("$dmap/*.male.map");
open Out,">$dOut/total.male.map";
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	print Out "group\t$lgID\n";
	open In,$map;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		print Out $_,"\n";
	}
	close In;
}
close Out;
@map=glob("$dmap/*.female.map");
open Out,">$dOut/total.female.map";
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	print Out "group\t$lgID\n";
	open In,$map;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		print Out $_,"\n";
	}
	close In;
}
close Out;

my @marker=glob("$dmap/*.correct.loc");
open Out,">$dOut/total.loc";
my $nloc=0;
my $nind=0;
my @out;
foreach my $marker (@marker) {
	open In,$marker;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my @info=split;
		next if (scalar @info < 3);
		if (/=/) {
			next;
		}else{
			if ($nind != scalar @info -3 && $nind != 0) {
				die "error $_";
			}else{
				$nind=scalar @info-3;
			}
			$nloc++;
			push @out,$_;
		}
	}
	close In;
}
print Out join("\n","nloc = $nloc","nind = $nind","popt = CP","name = Total",@out);
close Out;
open Out,">$dOut/total.trt";
print Out "ntrt=1\n";
print Out "nind=$nind\n";
print Out "miss=*\n";
print Out "T1\n";
for (my $i=0;$i<$nind;$i++) {
	print Out "1\n";
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
  -dmap	<file>	input file name
  -out	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}
