#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$dsh,$proc);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"out:s"=>\$outdir,
	"dsh:s"=>\$dsh,
	"proc:s"=>\$proc,
			) or &USAGE;
&USAGE unless ($fqlist and $outdir and $dsh);
$proc||=20;
mkdir $outdir if (!-d $outdir);
$outdir = ABSOLUTE_DIR($outdir);
mkdir $dsh if(!-d $dsh);
$dsh = ABSOLUTE_DIR($dsh);
open SH,">$dsh/step03.ustacks.sh";
open Out,">$outdir/ustacks.list";
open In,$fqlist;
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$fq)=split(/\s+/,$_);
	print SH "ustacks -f $fq -o $outdir/ -i $n --deleverage -M 6 -m 2 -p 8 --gapped\n";
	print Out "$sample\t$outdir/$sample\n";
	$n++;
}
close In;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=10G --CPU 8 --maxjob $proc $dsh/step03.ustacks.sh";
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
  -fqlist	<file>	input fqlist file
  -out	<dir>	split windows sh
  -dsh	<dir>	output work sh	
  -proc	<num>	max proc number for qsub
  -h         Help

USAGE
        print $usage;
        exit;
}
