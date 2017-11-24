#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dOut,$proc,$dShell,$ref,$dict,$split);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"ref:s"=>\$ref,
	"dict:s"=>\$dict,
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
$dict=ABSOLUTE_DIR($dict);
$split||=20;
open SH,">$dShell/07.haplotyper.sh";
open In,$bamlist;
open List,">$dOut/gvcf.list";
my %bam;
my $number=0;
my $nct=8;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	my $bai=$bam;
	$bai=~s/bam$/bai/;
	open Out,">$dOut/$sampleID.json\n";
	print Out "{\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.RefFasta\": \"$ref\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.Refdict\": \"$dict\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.Refindex\": \"$ref.fai\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.workdir\": \"$dOut\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.inputBAM\": \"$bam\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.sampleName\": \"$sampleID\",\n";
	print Out "\"HaplotypeCaller.haplotypeCaller.BAMindex\": \"$bai\"\n";
	print Out "}\n";
	close Out;
	print List join("\t",$sampleID,"$dOut/$sampleID.gvcfs"),"\n";
;	print SH "java -jar /mnt/ilustre/users/dna/.env//bin//cromwell-29.jar run $Bin/bin/HaplotypeCaller.wdl -i $dOut/$sampleID.json --metadata-output $dOut/\n";
}
close In;
close SH;
close List;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=30G --CPU 8 --maxjob $proc $dShell/07.haplotyper.sh";
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
  -dict	<file>	input dict file
  -proc <num>	number of process for qsub,defaulr
  -out	<dir>	output dir
  -dsh	<dir>	output shell dir
  -split	<num>	split files for speed default 20
  -h         Help

USAGE
        print $usage;
        exit;
}
