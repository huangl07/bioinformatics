#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$dOut,$dShell);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($fqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if(!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$dShell=ABSOLUTE_DIR($dShell);
open SH,">$dShell/step00.uniform.sh";
open Out,">$dOut/fq.uniform.list";
open In,$fqlist;
my %fq1;
my %fq2;
my %sample;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($sample,$fq1,$fq2)=split(/\s+/,$_);
	if (!-f $fq1 ||!-f $fq2) {
		die "$fq1|$fq2";
	}
	$sample{$sample}=1;
	push @{$fq1{$sample}},$fq1;
	push @{$fq2{$sample}},$fq2;
}
close In;
foreach my $sample (sort keys %sample) {
	my $fq1=join(" ",@{$fq1{$sample}});
	my $fq2=join(" ",@{$fq2{$sample}});
	print SH "cat $fq1 > $dOut/$sample.total.R1.fastq.gz && ";
	print SH "cat $fq2 > $dOut/$sample.total.R2.fastq.gz && ";
	print SH "$Bin/bin/fastq_length_uniform -1 $dOut/$sample.total.R1.fastq.gz -2 $dOut/$sample.total.R2.fastq.gz -l 120 -a $dOut/$sample\_R1.fastq.gz -b $dOut/$sample\_R2.fastq.gz && ";
	print SH "ln -s $dOut/$sample\_R1.fastq.gz $dOut/$sample.fastq.gz \n";
	print Out "$sample\t$dOut/$sample.fastq.gz\n";
}
close Out;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=3G --CPU 1 --Nodes 1 $dShell/step00.uniform.sh";
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
  -fqlist	<file>	input file name
  -out	<dir>	output dir
  -dsh	<dir>	workshell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
