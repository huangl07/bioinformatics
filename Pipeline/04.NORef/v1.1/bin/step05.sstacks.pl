#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ulist,$clist,$dOut,$dShell,$proc);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ulist:s"=>\$ulist,
	"clist:s"=>\$clist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
	"proc:s"=>\$proc,
			) or &USAGE;
&USAGE unless ($ulist and $clist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dShell if(!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open In,$clist;
open SH,">$dShell/step05.sstacks.sh";
open Out,">$dOut/sstacks.list";
my $cstack;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$cstacks)=split(/\s+/,$_);
	$cstack="$cstacks";
}
close In;
open In,$ulist;
my @ustacks;
print $cstack;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$ustacks)=split(/\s+/,$_);
	print SH "/mnt/ilustre/users/dna/.env/stacks-2.1/bin/sstacks -s $ustacks -c $cstack -p 15 -o $dOut/\n";
	print Out "$sample\t$dOut/$sample\n";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=100G --CPU 15 --Nodes 1 --Maxjob $proc $dShell/step05.sstacks.sh";
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
  -ulist	<file>	input ulist
  -clist	<file>	input clist
  -out	<dir>	output dir
  -dsh	<dir>	output worksh
  -h         Help

USAGE
        print $usage;
        exit;
}
