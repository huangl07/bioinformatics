#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamlist,$dOut,$dShell,$ref,$proc,$dict,$split);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"dict:s"=>\$dict,
	"ref:s"=>\$ref,
	"out:s"=>\$dOut,
	"proc:s"=>\$proc,
	"dsh:s"=>\$dShell,
	"split:s"=>\$split,
			) or &USAGE;
&USAGE unless ($bamlist and $dOut and $dShell);
$ref=ABSOLUTE_DIR($ref);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$bamlist=ABSOLUTE_DIR($bamlist);
$proc||=20;
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$split||=20;
#open In,$dict;
#my %len;
#while (<In>) {
	#chomp;
	#next if ($_ eq ""||/^$/);
	#my (undef,$sca,$length,undef)=split(/\s+/,$_);
	#$sca=~s/SN://g;
	#$length=~s/LD://g;
	#$len{$sca}=$length;
#}
#close In;
#my @sca=scalar keys %sca;
#my %bed;
#for (my $i=0;$i<@sca;$i++) {
	#my $hand=$i/$split;
	#if (!exists $bed{$hand}) {
		#open $bed{$hand},">$dOut/$hand.bed";
	#}
	#print {$bed{$hand}} join("\t",$sca[$i],0,$len{$sca}),"\n";
#}
#foreach my $head (sort keys %bed) {
	#close $bed{$hand};
#}


open SH,">$dShell/step05.call-variant.sh";
open In,$bamlist;
open Out,">$dOut/bam.list";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sampleID,$bam)=split(/\s+/,$_);
	print Out $bam,"\n";
}
close In;
close Out;
print SH "samtools mpileup  -t DP,AD -uf $ref -b $dOut/bam.list|bcftools call -mv --format-fields GQ,GP --output-type z > $dOut/var.raw.vcf.gz ";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=30G --CPU 1 --maxjob $proc $dShell/step05.call-variant.sh";
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
	perl $Script -bam -dict -ref -out -dsh 

Usage:
  Options:
  -bam	<file>	input bamlist file
  -ref	<file>	input reference file
  -out	<dir>	output dir
  -proc <num>	number of process for qsub,default
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
