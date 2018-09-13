#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gvcflist,$dOut,$proc,$dShell,$ref);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gvcf:s"=>\$gvcflist,
	"ref:s"=>\$ref,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($gvcflist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$ref=ABSOLUTE_DIR($ref);
my $sentieon="/mnt/ilustre/users/dna/.env/sentieon/sentieon/sentieon-genomics-201711.05/bin/sentieon";
open SH,">$dShell/07.gvcf-typing.sh";
print SH "$sentieon driver -r $ref --algo GVCFtyper --var_type both ";
open In,$gvcflist;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$gvcf)=split(/\s+/,$_);
	if (!-f $gvcf|| !-f "$gvcf.idx") {
		die "check $gvcf!";
	}
	print SH "-v $dOut/$sampleID.g.vcf ";
}
print SH "$dOut/pop.variant.vcf\n";
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=30G --CPU 1 $dShell/07.gvcf-typing.sh";
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
  -gvcf	<file>	input gvcflist file
  -ref	<file>	input reference file
  -out	<dir>	output dir
  -dsh	<dir>	output shell dir

  -h         Help

USAGE
        print $usage;
        exit;
}
