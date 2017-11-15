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
$fqlist=ABSOLUTE_DIR($fqlist);
open In,$fqlist;
open SH,">$dShell/step05.clean-generic.sh";
open Out,">$dOut/fastq.list";
my %lane;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$fq1,$fq2)=split(/\s+/,$_);
	print SH "java -Xmx16G -jar /mnt/ilustre/users/dna/.env//bin/trimmomatic-0.36.jar  PE -threads 8 -phred33 -trimlog $dOut/$sample.ada.log ";
	print SH "$fq1 $fq2 $dOut/$sample.clean.1.fastq.gz $dOut/$sample.unpair.1.fq.gz $dOut/$sample.clean.2.fastq.gz  $dOut/$sample.unpair.2.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n";
	print Out "$sample $dOut/$sample.clean.1.fastq.gz $dOut/$sample.clean.2.fastq.gz\n";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dShell/step05.clean-generic.sh --Resource mem=3G";
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
	-fqlist	<file>	input fqlist file 
	-out	<file>	output dir
	-dsh	<file>	output work shell
USAGE
	print $usage;
	exit;
}
