#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$ref,$gff,$RAD,$step,$stop,$SV,$CNV,$realign);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"anno:s"=>\$anno,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $out and $anno);
mkdir $out if (!-d $out);
$vcf=ABSOLUTE_DIR($vcf);
$anno=ABSOLUE_DIR($anno);
my $dsh="$outdir/work_sh";
mkdir $dsh if (!-d $dsh);
open LOG,">$outdir/work_sh/BSA.$BEGIN_TIME.log";
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "index calc\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastq-qc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "sliding window\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastq-qc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "region abstract\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastq-qc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
close LOG;
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
	eg:
	perl $Script -fqlist -out -ref -gff

Usage:

	-fqlist	<file>	input file name
	-outdir	<dir>	output dir
	-ref	<file>	reference file
	-gff	<file>	gff file
	-sv	sv calling default off 
	-cnv	cnv calling default off 
	-RAD	RRL calling default off 
	
	-step	pipeline control
          01 fastq qc
          02 reference prepair
          03 bwa mapping
          04 bam sort
          05 bam mkdup
          06 map-stat
          07 haplotype
          08 gvcf typing
          09 vcf-filter
          10 annovar
          11 sv call
          12 cnv call
          13 variant stat
          14 report
	-stop	pipeline control

  -h         Help

USAGE
        print $usage;
        exit;
}
