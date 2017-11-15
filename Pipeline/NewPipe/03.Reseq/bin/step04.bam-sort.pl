#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$outdir,$proc,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"outdir:s"=>\$outdir,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dsh
			) or &USAGE;
&USAGE unless ($bamlist and $outdir and $dsh);
mkdir $outdir if (!-d $outdir);
$outdir=ABSOLUTE_DIR($outdir);
$bamlist=ABSOLUTE_DIR($bamlist);
$proc||=20;
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open SH,">$dsh/04.bam-sort.sh";
open Out,">$outdir/bam.sort.list";
open In,$bamlist;
my %bam;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,@bam)=split(/\s+/,$_);
	foreach my $bam (@bam) {
		if (!-f $bam) {
			die "check $bam!";
		}
	}
	my $bam=join(" I=",@bam);
	print Out $sampleID,"\t$dOut/$sampleID.sort.bam\n";
	print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/picard.jar MergeSamFiles I=$bam O=$dOut/$sampleID.sort.bam SORT_ORDER=coordinate TMP_DIR=$dOut/merge MAX_RECORDS_IN_RAM=50000000 VALIDATION_STRINGENCY=LENIENT&& ";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=20G --CPU 1 --maxjob $proc  $dsh/04.bam-sort.sh";
`$job`;
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
	perl $Script -bam -out -dsh

Usage:
  Options:
  -bam	<file>	input bamlist file
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
