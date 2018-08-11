#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$dsh,$GBS,$proc,$method);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"method:s"=>\$method,
	"out:s"=>\$outdir,
	"dsh:s"=>\$dsh,
	"proc:s"=>\$proc,
			) or &USAGE;
&USAGE unless ($fqlist and $outdir and $dsh);
$proc||=20;
mkdir $outdir if (!-d $outdir);
mkdir $dsh if(!-d $dsh);
$outdir=ABSOLUTE_DIR($outdir);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step02.uniform.sh";
open Out,">$outdir/fq.list";
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
	print SH "cat $fq1 > $outdir/$sample.total.R1.fastq.gz && ";
	print SH "cat $fq2 > $outdir/$sample.total.R2.fastq.gz && ";
	print SH "$Bin/bin/ngsuniform -1 $outdir/$sample.total.R1.fastq.gz -2 $outdir/$sample.total.R2.fastq.gz -l 120 -a $outdir/$sample\_R1.fastq.gz -b $outdir/$sample\_R2.fastq.gz && ";
	if ($method eq "GBS") {
		print SH "cat $outdir/$sample\_R1.fastq.gz $outdir/$sample\_R2.fastq.gz > $outdir/$sample.fastq.gz \n";
	}else{
		print SH "ln -s $outdir/$sample\_R1.fastq.gz $outdir/$sample.fastq.gz \n";
	}
	print Out "$sample\t$outdir/$sample.fastq.gz\n";
}
close Out;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue dna --Resource mem=3G --CPU 1 --maxjob $proc $dsh/step02.uniform.sh";
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
  -method	<str>	GBS or RAD default RAD
  -dsh	<dir>	workshell dir
  -proc	<num>	max proc number for qsub
  -h         Help

USAGE
        print $usage;
        exit;
}
