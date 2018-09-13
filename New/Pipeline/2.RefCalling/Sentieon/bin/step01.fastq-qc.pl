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
mkdir "$outdir/stat/" if (!-d "$outdir/stat/");
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
my %sample;
my %fq1;
my %fq2;
open In,$fqlist;
open SH,">$dsh/01.fastq-qc.sh";
my %nsample;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$fq1,$fq2,$type)=split(/\s+/,$_);
	if (!defined $type) {
		die "error fqlist file! add library type";
	}
	if (!-f $fq1 || !-f $fq2) {
		die "!exists $fq1 | $fq2";
	}
	$sample{$id}=$type;
	push @{$fq1{$id}},$fq1;
	push @{$fq2{$id}},$fq2;
}
close In;
open OUT,">$outdir/fq.list";
foreach my $id (sort keys %sample) {
	my $fq1=join(" ",@{$fq1{$id}});
	my $fq2=join(" ",@{$fq2{$id}});
	print OUT "$id\t$outdir/$id.total.R1.fastq.gz\t$outdir/$id.total.R2.fastq.gz\t$sample{$id}\n";
	print SH "cat $fq1 > $outdir/$id.total.R1.fastq.gz && ";
	print SH "cat $fq2 > $outdir/$id.total.R2.fastq.gz && ";
	print SH "ngsqc -1 $outdir/$id.total.R1.fastq.gz -2 $outdir/$id.total.R2.fastq.gz -o $outdir/stat -k $id  && ";
	print SH "Rscript $Bin/bin/ngsqc.R --base $outdir/stat/$id.atgc --qual $outdir/stat/$id.qual --out $outdir/fig/$id\n";
}
close OUT;
close SH;
my $job="perl  /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10G --CPU 1 --maxjob $proc $dsh/01.fastq-qc.sh";
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
Contact: caixia.tian

Usage:

	-fqlist	<file>	input fastq list	
	-outdir	<dir>	output result dir
	-dsh	<dir>	output work dir
	-proc	<num>	max proc number for qsub

USAGE
	print $usage;
	exit;
}
