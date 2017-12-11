#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fqlist,$config,$dOut,$dShell);
GetOptions(
				"help|?" =>\&USAGE,
				"fqlist:s"=>\$fqlist,,
				"out:s"=>\$dOut,
				"dsh:s"=>\$dShell,
				) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
mkdir "$dOut/fig" if (!-d "$dOut/fig");
$dOut=ABSOLUTE_DIR($dOut);
$fqlist=ABSOLUTE_DIR($fqlist);
open In,$fqlist;
open SH,">$dShell/step02.rawdata-qc.sh";
open Out,">$dOut/stat.list";
mkdir "$dOut/fig" if (!-d "$dOut/fig");
my %lane;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2,$e1,$e2)=split(/\s+/,$_);
	print SH "$Bin/bin/ngsqc -1 $fq1 -2 $fq2 -o $dOut -k $sample -e $e1 -p $e2 && ";
	print SH "Rscript $Bin/bin/ngsqc.r --base  $dOut/$sample.atgc  --qual $dOut/$sample.qual  --key $sample --od  $dOut/fig\n";
	print Out "$sample $dOut/$sample.stat\n";
}
close In;
close SH;
close Out;
my $job="perl $Bin/bin/qsub-sge.pl $dShell/step02.rawdata-qc.sh --Resource mem=3G --CPU 1  --Nodes 1 --Maxjob 19 ";
`$job`;

#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-fqlist	<file>	input file 
	-out	<dir>	output dir
	-dsh	<dir>	output dsh
USAGE
	print $usage;
	exit;
}
