#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dOut,$dShell,$ref,$dict,$proc,$split);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"ref:s"=>\$ref,
	"dict:s"=>\$dict,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell,
	"split:s"=>\$split,
			) or &USAGE;
&USAGE unless ($bamlist and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$bamlist=ABSOLUTE_DIR($bamlist);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$split||=20;
open In,$dict;
my %chr;
my $nsplit=1;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || !/^\@SQ/);
	my ($SQ,$SN,$LN,undef)=split(/\s+/,$_);
	my $id=(split(/\:/,$SN))[-1];
	$nsplit++;
	my $n=$nsplit % $split;
	$chr{$n}.=" -L $id ";
}
close In;
open SH,">$dShell/step05.Haplotyper.sh";
open In,$bamlist;
open Out,">$dOut/gvcf.list";
my %bam;
my $number=0;
my $nct=8;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	foreach my $chr (sort keys %chr) {
		print Out join("\t",$sampleID,$chr,"$dOut/$sampleID.$chr.g.vcf"),"\n";
		print SH "java -Djava.io.tmpdir=$dOut/tmp/ -Xmx30G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -nct $nct -R $ref -I $bam --genotyping_mode DISCOVERY -o $dOut/$sampleID.$chr.g.vcf -log $dOut/$sampleID.$chr.HaplotypeCaller.log $chr{$chr} --emitRefConfidence GVCF -stand_call_conf 30 -variant_index_type LINEAR -variant_index_parameter 128000 -filterNoBases -filterMBQ -filterRNC\n";
	}
}
close In;
close SH;
close Out;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=30G --CPU 8 --maxjob $proc  $dShell/step05.Haplotyper.sh";
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
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	output shell dir
  -split	<num>	split files for speed default 20
  -h         Help

USAGE
        print $usage;
        exit;
}
