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
	"out:s"=>\$outdir,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dsh
			) or &USAGE;
&USAGE unless ($bamlist and $outdir and $dsh);
mkdir $outdir if (!-d $outdir);
$outdir=ABSOLUTE_DIR($outdir);
$bamlist=ABSOLUTE_DIR($bamlist);
$proc||=20;
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/04.bam-sort.sh";
open Out,">$outdir/bam.list";
open In,$bamlist;
my %bam;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$type,@bam)=split(/\s+/,$_);
	foreach my $bam (@bam) {
		if (!-f $bam) {
			die "check $bam!";
		}
	}
	my $bam=join(" ",@bam);
	print Out $sampleID,"\t$type\t$outdir/$sampleID.sort.bam\n";
	print SH "samtools merge -f -p -@ 8 --output-fmt BAM $outdir/$sampleID.merged.bam $bam && ";
	print SH "samtools sort -o $outdir/$sampleID.sort.bam --output-fmt BAM -@ 8 $outdir/$sampleID.merged.bam &&";
	print SH "samtools index $outdir/$sampleID.sort.bam \n";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=20G --CPU 8 --maxjob $proc $dsh/04.bam-sort.sh";
print $job;
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
