#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$gvcflist,$dOut,$dShell,$ref,$dict);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gvcf:s"=>\$gvcflist,
	"ref:s"=>\$ref,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell
			) or &USAGE;
&USAGE unless ($gvcflist and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$gvcflist=ABSOLUTE_DIR($gvcflist);
$ref=ABSOLUTE_DIR($ref);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
open SH,">$dShell/08.gvcf-typing.sh";
open In,$gvcflist;
my $vcf;
my $number=0;
my $nct=8;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$gvcf)=split(/\s+/,$_);
	$vcf .=" -V $gvcf ";
}
close In;
print SH "java -Djava.io.tmpdir=$dOut/tmp/ -Xmx35G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T GenotypeGVCFs $vcf -nt 8 -o $dOut/pop.noid.vcf -R $ref -log $dOut/pop.variant.log && ";
print SH "bcftools annotate --set-id +\'\%CHROM\\_\%POS\' $dOut/pop.noid.vcf -o $dOut/pop.variant.vcf\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=35G --CPU 8 --maxjob $proc $dShell/08.gvcf-typing.sh";
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
  -gvcf	<file>	input bamlist file
  -ref	<file>	input reference file
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
