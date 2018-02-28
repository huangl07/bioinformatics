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
my ($fqlist,$dOut,$dShell);
GetOptions(
				"help|?" =>\&USAGE,
				"fqlist:s"=>\$fqlist,,
				"out:s"=>\$dOut,
				"dsh:s"=>\$dShell,
				) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
`mkdir $dOut/qc/ `;
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);

open In,$fqlist;
open SH,">$dShell/step03.fastp-trim.sh";
mkdir "$dOut/fig" if (!-d "$dOut/fig");
open Out,">$dOut/fastq.list";
open Stat,">$dOut/stat.list";
open Json,">$dOut/json.list";
my %lane;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2,$e1,$e2)=split(/\s+/,$_);
	print SH "/mnt/ilustre/users/dna/.env//bin/fastp  -i $fq1 -I $fq2 -o $dOut/$sample.clean.1.fastq.gz -O $dOut/$sample.clean.2.fastq.gz ";
	print SH " -q 20 -l 36 -5 20 -3 3 -W 4 -M 20 -n 10  -j $dOut/qc/$sample.json -h $dOut/qc/$sample.html -z 6 -w 8 && ";
	print SH " perl /mnt/ilustre/users/caixia.tian/perl/fastp/bin/bin/fastp.pl -i $dOut/qc/$sample.json -o $dOut/qc/$sample && Rscript /mnt/ilustre/users/caixia.tian/perl/fastp/bin/bin/ngsqc.r --base  $dOut/qc/$sample.raw.atgcn  --qual $dOut/qc/$sample.raw.qual  --key $sample --od  $dOut/fig && Rscript /mnt/ilustre/users/caixia.tian/perl/fastp/bin/bin/ngsqc.r --base $dOut/qc/$sample.clean.atgcn --qual $dOut/qc/$sample.clean.qual --key $sample --od $dOut/fig \n";
	print Out "$sample $dOut/$sample.clean.1.fastq.gz $dOut/$sample.clean.2.fastq.gz\n";
	print Stat "$sample\t$dOut/qc/$sample.stat\n";
	print Json "$sample\t$dOut/qc/$sample.json\n";
}
close In;
close SH;
close Out;
close Stat;
close Json;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dShell/step03.fastp-trim.sh --Resource mem=12G --CPU 8  --Maxjob 39";
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
	-dsh	<dir>	output work sh
USAGE
	print $usage;
	exit;
}
