#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bam,$out,$proc,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bam,
	"out:s"=>\$out,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($bam and $out and $dsh);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$bam=ABSOLUTE_DIR($bam);
$proc||=20;
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step02.bam-sort.sh";
open Out,">$out/bam.sort.list";
open In,$bam;
my %bam;
my $number=0;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,@bam)=split(/\s+/,$_);
	foreach my $bam (@bam) {
		if (!-f $bam) {
			die "check $bam!";
		}
	}
	my $bam=join(" ",@bam);
	$number++;
	print Out $sampleID,"\t$out/$sampleID.sort.bam\n";
	print SH "samtools merge -f -p -@ 8 --output-fmt BAM $out/$sampleID.merged.bam $bam && ";
	print SH "samtools sort -o $out/$sampleID.sort.bam --output-fmt BAM -@ 8 $out/$sampleID.merged.bam &&";
	print SH "samtools index $out/$sampleID.sort.bam ";
	print SH "\n";
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=20G --CPU 1 --maxjob $proc  $dsh/step02.bam-sort.sh";
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
  -bam	<file>	input bamlist
  -out	<dir>	ouput dir
  -proc <num>   number of process for qsub,default 20
  -dsh	<dir>	output work_sh shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
