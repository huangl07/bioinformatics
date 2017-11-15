#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ulist,$clist,$dOut,$dShell,$slist);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ulist:s"=>\$ulist,
	"group:s"=>\$group,
	"clist:s"=>\$clist,
	"slist:s"=>\$slist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($ulist and $clist and $slist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if(!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$dShell=ABSOLUTE_DIR($dShell);
mkdir "$dOut/stacks" if (!-d "$dOut/stacks");

open In,$ulist;
my $Austacks;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$ustacks)=split(/\s+/,$_);
	`ln -s $ustacks* $dOut/stacks`;
	my $ustack=basename($ustacks);
	$Austacks.=" -s $dOut/$ustack ";
}
close In;
open In,$clist;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$cstacks)=split(/\s+/,$_);
	`ln -s $cstacks* $dOut/stacks`;
}
close In;
open In,$slist;
open SH,">$dShell/step06.correct2.sh";
open List,">$dOut/sstacks.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$sstacks)=split(/\s+/,$_);
	`ln -s $sstacks* $dOut/stacks`;
	print SH "sstacks -b 1 -c $dOut/batch_1 -s $dOut/$sample -o $dOut -p 8 --gapped\n";
	print List "$sample\t$dOut/$sample\n";
}
close SH;
close List;
close In;
open SH,">$dShell/step06.correct1.sh";
print SH "rxstacks -P $dOut/stacks/ --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.1 --lnl_lim -10.0 -t 16 -o $dOut/ && ";
print SH "cstacks $Austacks/ -n 4 -p 16 --gapped -o $dOut\n";
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=80G --CPU 16 --Nodes 1 $dShell/step06.correct1.sh";
print  "$job\n";
`$job`;
print "$job\tdone!\n";
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=80G --CPU 16 --Nodes 1 $dShell/step06.correct2.sh";
print  "$job\n";
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
	-ulist	<file>	ustacks list
	-clist	<file>	cstacks	list
	-slist	<file>	sstacks list
	-group	<file>	group list for cstacks
	-out	<dir>	output dir
	-dsh	<dir>	work shell dir
	-h         Help

USAGE
        print $usage;
        exit;
}
