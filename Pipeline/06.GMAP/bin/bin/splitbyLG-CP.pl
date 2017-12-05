#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fLG,$fKey,$type);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$fLG,
				"d:s"=>\$dOut,
				"t:s"=>\$type,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $fLG);
open In,$fIn;
my %info;
my $head;
my $nind;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		$head=$_;
		next;
	}
	my ($id,$info)=split(/\s+/,$_,2);
	$info{$id}=$info;
	my @ind=split(/\t/,$info);
	$nind=scalar @ind;
}
close In;
open In,$fLG;
open List,">$dOut/pri.marker.list";
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,$marker)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	open Out,">$dOut/$id.marker";
	print List "$id\t$dOut/$id.marker\n";
	my @marker=split(/\s+/,$marker);
	my @out;
	my $nloc=scalar @marker;
	foreach my $m (@marker) {
		push @out,join("\t",$m,$info{$m});
	}
	print Out $head,"\n";
	print Out join("\n",@out),"\n";
	close Out;
}
close In;
close List;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: huangl <long.huang\@majorbio.com> 

Usage:  将genotype文件按连锁群分割
  Options:
  -help			USAGE,
  -i	genotype file， forced
  -l	linkage lg file
  -o	output dir
  
   
USAGE
	print $usage;
	exit;
}

