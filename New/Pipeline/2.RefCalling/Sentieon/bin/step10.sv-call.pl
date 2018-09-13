#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$bamlist,$dOut,$dShell,$ref,$con);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"ref:s"=>\$ref,
	"con:s"=>\$con,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell
			) or &USAGE;
&USAGE unless ($bamlist and $ref and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$bamlist=ABSOLUTE_DIR($bamlist);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
my $sentieon="/mnt/ilustre/users/dna/.env/sentieon/sentieon/sentieon-genomics-201711.05/bin/sentieon";
open SH1,">$dShell/10.sv-calling1.sh";
open SH2,">$dShell/10.sv-calling2.sh";
print SH2 "$sentieon driver -t 8 -r $ref --algo SVSolver ";
open In,$bamlist;
open List,">$dOut/sv.list";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	if (!-f $bam) {
		die "check $bam!";
	}
	print SH1 "$sentieon driver -t 8 -r $ref -i $bam --algo DNAscope --var_type bnd $dOut/$sampleID.tmp.vcf \n";
	print SH2 "-v $dOut/$sampleID.tmp.vcf ";
	#print List "$sampleID $dOut/pop.sv.vcf\n";
}
print List "$sampleID $dOut/sv.anno.primary.vcf\n";
print SH2 "$dOut/pop.sv.vcf && java -Xmx50G -jar /mnt/ilustre/users/dna/.env//bin//snpEff.jar -v ref -csvStats $dOut/sv.anno.csv -c $con $dOut/pop.sv.vcf > $dOut/sv.anno.primary.vcf && perl $Bin/bin/sv.stat.pl -i $dOut/sv.anno.csv -o $dOut/ \n";
close In;
close SH1;
close SH2;
close List;

my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=20G --CPU 1  --maxjob $proc $dShell/10.sv-calling1.sh";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=50G --CPU 8  --maxjob $dShell/10.sv-calling2.sh";
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
  -ref	<file>	input refrence file
  -con	<file>	input snpEff.config
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default 20
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
