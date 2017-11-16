#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$vcf,$dOut,$dShell,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"proc:s"=>\$proc,
	"vcf:s"=>\$vcf,
	"ref:s"=>\$ref,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($vcf and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$vcf=ABSOLUTE_DIR($vcf);
$ref=ABSOLUTE_DIR($ref);
open SH,">$dShell/09.vcf-filter.sh";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V $vcf -selectType SNP  -o $dOut/pop.snp.vcf -log $dOut/pop.selectSNP.log && ";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V $dOut/pop.snp.vcf --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SQR > 3.0\" --filterName Failer -o $dOut/pop.snp.filter.vcf -log $dOut/pop.snp.filter.log&&";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V $dOut/pop.snp.filter.vcf -o $dOut/pop.snp.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered --excludeNonVariants -log $dOut/pop.filterSNP.log \n";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T SelectVariants  -R $ref -V $vcf -selectType INDEL  -o $dOut/pop.indel.vcf -log $dOut/pop.selectINDEL.log && ";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V $dOut/pop.indel.vcf  --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0\" --filterName Failer -o $dOut/pop.indel.filter.vcf -log $dOut/pop.indel.filter.log&& ";
print SH "java -Xmx20G -jar /mnt/ilustre/users/dna/.env/bin/GenomeAnalysisTK.jar -T SelectVariants  -R $ref -V $dOut/pop.indel.filter.vcf -o $dOut/pop.indel.filter.recode.vcf --setFilteredGtToNocall --excludeFiltered --excludeNonVariants -log $dOut/pop.selectINDEL.log \n";
close SH;
open List,">$dOut/vcf.filter.list";
print List "snp\t$dOut/pop.snp.filter.recode.vcf\n";
print List "indel\t$dOut/pop.indel.filter.recode.vcf\n";
close List;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=20G --CPU 1  --maxjob $proc $dShell/09.vcf-filter.sh";
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
	fq thanslate to ref format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input genome vcf file
  -ref	<file>	input reference file
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	work sh
  -h         Help

USAGE
        print $usage;
        exit;
}
