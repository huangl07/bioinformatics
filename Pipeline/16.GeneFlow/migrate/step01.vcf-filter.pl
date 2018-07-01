#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$out,$dsh,$maf,$mis,$dep,$gro);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
    "group:s"=>\$group,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh );
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
$group=ABSOLUTE_DIR($group);
$mis||=0.3;
$maf||=0.05;
$dep||=2;
$mis=1-$mis;
open Out,">$out/pop.parmfile";
my $parm=<<"PARM";
################################################################################
# Parmfile for Migrate debug 3.6.6-Dec15/15 [do not remove these first TWO lines]
menu=no
nmlength=30
datatype=SequenceData
ttratio=2.000000 
freqs-from-data=YES
seqerror-rate=0.000000
rates=1: 1.000000 
prob-rates=1: 1.000000 
autocorrelation=NO
weights=NO
interleaved=NO
fast-likelihood=NO
inheritance-scalars={1.00000000000000000000}
population-relabel={1}
usertree=AUTOMATIC
infile=pop.migrate
title=pop analysis
progress=YES
logfile=YES:logfile
print-data=NO
outfile=pop.migrate.out
pdf-outfile=pop.migrate.out.pdf
use-M=YES
plot=NO
mathfile=pop.mathfile
profile=ALL:PRECISE
print-tree=BEST:pop.best.tree
write-summary=YES
aic-modeltest=NO
mig-histogram=YES:0.001000:mighisfile
theta=FST
migration=FST
mutation=ESTIMATE
fst-type=THETA
custom-migration={**}
geo=NO
bayes-update=YES
bayes-updatefreq=0.500000
bayes-posteriorbins=1500 1500
bayes-posteriormaxtype=TOTAL
bayes-file=NO
bayes-allfile=NO
bayes-proposals= THETA METROPOLIS-HASTINGS Sampler
bayes-proposals= MIG METROPOLIS-HASTINGS Sampler
bayes-priors= THETA UNIFORMPRIOR: 0.000000 0.100000 0.010000 
bayes-priors= MIG UNIFORMPRIOR: 0.000000 1000.000000 100.000000 
long-chains=1
long-inc=100
long-sample=5000
burn-in=10000  
auto-tune=YES:0.440000
heating=YES:1:{1.000000,1.500000,3.000000,1000000.000000}
heated-swap=YES
moving-steps=NO
long-chain-epsilon=INFINITY
gelman-convergence=No
replicate=NO
resistance=0.000001
PARM
print Out $parm;
close Out;
open SH,">$dsh/migrate.sh";
print SH "vcftools --remove-filtered-all --remove-indels --minDP $dep  --max-missing $mis --vcf $vcf --recode  --out $out/pop --maf $maf &&";
print SH "perl $Bin/vcf2migrate.pl -vcf $out/pop.recode.vcf -group $group -out $out/pop.migrate && ";
print SH "cd $out/ && mpirun -np 32 /mnt/ilustre/users/dna/.env/bin/migrate-n-mpi $out/pop.parmfile";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --CPU 32 --Resource mem=20G $dsh/migrate.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	filter snp and retrive fa 
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf files
  -group <file> input group list;split by \\t
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
