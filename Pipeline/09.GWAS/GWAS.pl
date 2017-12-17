#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro,$trt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"trt:s"=>\$trt,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
	"block:s"=>\$block,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh and $block);
my $vcf=ABSOLUTE_DIR($vcf);
my $trt=ABSOLUTE_DIR($trt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$mis||=0.3;
$maf||=0.05;
$dep||=2;
$mis=1-$mis;

open SH,">$dsh/gwas1.sh";
print SH "vcftools --vcf $vcf --out $out/myData --min-meanDP $dep --max-missing $mis --maf $maf --recode && ";
print SH "ln -s $out/myData.recode.vcf $out/myData.vcf && ";
print SH "ln -s $trt $out/myData.txt && cd $out/"
print SH "blink --file MyData --vcf  --gwas --trait 0 --parallel 8 \n";
close SH;
open SH,">$dsh/gwas2.sh";
open In,$trt;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@trt)=split(/\s+/,$_);
	for (my $i=0;$i<@trt;$i++) {
		print SH "Rscript $bin/manhattan.R --infile $out/$trt[$i]_GWAS_Result.txt --outfile  $out/$trt[$i]\n";
		print SH "perl get-region.pl -i  $out/$trt[$i]_GWAS_Result.txt -b $block -outfile  $out/$trt[$i].region\n";
	}
}
close In;
close SH;


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
  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -gro	<str>	group list
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2
  -block	<file>	block file

  -h         Help

USAGE
        print $usage;
        exit;
}
