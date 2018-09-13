#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dOut,$proc,$dShell,$ref,$dict,$split,$chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"ref:s"=>\$ref,
	"proc:s"=>\$proc,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($bamlist and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$bamlist=ABSOLUTE_DIR($bamlist);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$ref=ABSOLUTE_DIR($ref);
my $sentieon="/mnt/ilustre/users/dna/.env/sentieon/sentieon/sentieon-genomics-201711.05/bin/sentieon";
open SH,">$dShell/06.haplotyper.sh";
open In,$bamlist;
open List,">$dOut/gvcf.list";
my %bam;
my $number=0;
my $nct=8;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	if (!-f $bam) {
		die "check $bam!";
	}
	print SH "$sentieon driver -t 8 -r $ref -i $bam --algo Haplotyper $dOut/$sampleID.g.vcf\n";
	print List join("\t",$sampleID,"$dOut/$sampleID.g.vcf"),"\n";
}
close In;
close SH;
close List;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=30G --CPU 1 --maxjob $proc $dShell/06.haplotyper.sh";
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
	perl $Script -i -o -k -c

Usage:
  Options:
  -bam	<file>	input bamlist file
  -ref	<file>	input reference file
  -proc <num>	number of process for qsub,default 20
  -out	<dir>	output dir
  -dsh	<dir>	output shell dir

  -h         Help

USAGE
        print $usage;
        exit;
}
