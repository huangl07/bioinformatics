#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fqlist,$outdir,$dsh,$proc);
GetOptions(
				"help|?" =>\&USAGE,
				"fqlist:s"=>\$fqlist,
				"outdir:s"=>\$outdir,
				"dsh:s"=>\$dsh,
				"proc:s"=>\$proc,
				) or &USAGE;
&USAGE unless ($fqlist and $outdir and $dsh);
$proc||=20;
mkdir $outdir if (!-d $outdir);
$outdir=ABSOLUTE_DIR($outdir);
mkdir "$outdir/fig/" if (!-d "$outdir/fig");
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open In,$fqlist;
open SH,">$dsh/01.fastq-qc.sh";
my %nsample;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$fq1,$fq2)=split(/\s+/,$_);
	if (!-f $fq1 || !-f $fq2) {
		die "!exists $fq1 | $fq2";
	}
	$nsample{$id}++;
	print SH "ngsqc -1 $fq1 -2 $fq2 -o $outdir -k $id:$nsample{$id}  && ";
	print SH "Rscript $Bin/bin/ngsqc.R --base $outdir/$id-$nsample{$id}.atgc --qual $outdir/$id-$nsample{$id}.qual --out $outdir/fig/$id:$nsample{$id}\n";
}
close In;
close SH;
my $job="perl $Bin/../tools/qsub-sge.pl --Resource mem=3G --CPU 1 --maxjob $proc $dsh/01.fastq-qc.sh";
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

	-fqlist	<file>	input fastq list	
	-outdir	<dir>	output result dir
	-dsh	<dir>	output work dir
	-proc	<num>	max proc number for qsub

USAGE
	print $usage;
	exit;
}
