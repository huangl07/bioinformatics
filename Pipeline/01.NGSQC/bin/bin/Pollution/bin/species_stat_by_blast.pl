#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
	open In,$fIn;
	my %stat;
	my $sum=0;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my @info=split("\t",$_);
		my @ids=split(" ",$info[15]);
		my @id;
		for (my $i=1;$i<@ids;$i++) {
			next if ($ids[$i] =~ /sp./);
			next if ($ids[$i] =~ /PREDICTED:/);
			next if (length($ids[$i]) == 1);
			push @id,$ids[$i];
			last if (scalar @id ==2);
		}
		my $sample = join(" ",@id);
		$stat{$sample}++;
		$sum++;
	}
	close In;
	open Out,">$fOut";
	print Out "#Species\tNum\tPercentage\n";
	foreach my $out (sort{$stat{$b}<=>$stat{$a}} keys %stat) {
		print Out "$out\t$stat{$out}\t",$stat{$out},"\n";
	}
	close Out;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
