#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ulist,$clist,$dOut,$dShell,$sample);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ulist:s"=>\$ulist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
	"sample:s"=>\$sample,
			) or &USAGE;
&USAGE unless ($ulist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dShell if(!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
my %group;
if ($sample) {
	open In,$sample;
	open Out,">$dOut/sample.list";
	while (<In>) {
		chomp;
		next if ($_ eq "" ||/^$/);
		$group{$_}=1;
		print $_,"\n";
	}
	close In;
	close Out;
}else{
	open In,$ulist;
	my %dep;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($sample,$ustacks)=split(/\s+/,$_);
		my $total++;
		my $ndp++;
		open Slist,"gunzip -c $ustacks.tags.tsv.gz|";
		while (<Slist>) {
			chomp;
			next if ($_ eq ""||/^$/);
			if (/consensus/) {
				$total++;
			}
			if (/primary/ || /secondary/) {
				$ndp++;
			}
		}
		close Slist;
		$dep{$sample}=$ndp/$total;
	}
	close In;
	open Out,">$dOut/sample.list";
	my $n=0;
	foreach my $sam (sort {$dep{$a}<=> $dep{$b}}keys %dep) {
		if (scalar keys %dep > 10){
			next if ($n % 5 !=0);
		}
		$n++;
		$group{$sam}=1;
		print Out $sam,"\n";
	}
	close Out;
}
open In,$ulist;
open SH,">$dShell/step07.rcstacks.sh";
my $ustack;
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$ustacks)=split(/\s+/,$_);
	next if ($sample && !exists $group{$sample});
	if ($n == 0) {
		print SH "cstacks -b 1 -s $ustacks -n 4 -p 32 --gapped -o $dOut 2>$dOut/cstacks.$sample.log ";
	}else{
		print SH "&& cstacks -b 1 -s $ustacks -n 4 -p 32 --gapped -o $dOut --catalog $dOut/batch_1  2>$dOut/cstacks.$sample.log ";
	}
	$n++;
}
close In;
open Out,">$dOut/cstacks.list";
print Out "cstacks\t$dOut/batch_1\n";
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=20G --CPU 32 --Nodes 1 $dShell/step07.rcstacks.sh";
print "$job\n";
`$job`;
print "$job\tdone!\n";

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
	-ulist	<file>	input ustacks list file
	-out	<dir>	output dir
	-dsh	<dir>	output workshell dir
	-sample	<file>	sample list if no use all 
	-h         Help

USAGE
        print $usage;
        exit;
}
