#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($vcf,$out,$dsh,$pid,$bid,$qtl,$pop);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$vcf,
				"out:s"=>\$out,
				"dsh:s"=>\$dsh,
				"pid:s"=>\$pid,
				"bid:s"=>\$bid,
				"pop:s"=>\$pop,
				"qtl"=>\$qtl
				) or &USAGE;
&USAGE unless ($vcf and $out and $dsh);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$pop||="F2";
open SH,">$dsh/01.index-calc.sh";
open SH2,">$dsh/02.sliding-win.sh";
if ($qtl) {
	print SH "perl $Bin/bin/qtlseq.pl -vcf $vcf -out $out/index.calc -bid $bid -popt $pop ";
	if ($pid) {
		print SH "-pid $pid ";
	}
	print SH "Rscript --infile $out/index.calc --outfile $out/sliding.bp.result --col 1,2,14,15,16 --win 2000000 --col 10000 --method bp\n"
}else{
	print SH "per $Bin/bin/mutmap.pl -vcf $vcf -out $out/index.calc  -pid -bid $bid";
	if ($pid) {
		print SH "-pid $pid ";
	}
	print SH "Rscript --infile $out/index.calc --outfile $out/sliding.bp.result --col 1,2,10 --win 2000000 --col 10000 --method bp\n"
}
close SH;
close SH2;
my $job="perl $Bin/../tools/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/01.index-calc.sh";
`$job`;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
				"vcf:s"=>\$vcf,
				"out:s"=>\$out,
				"dsh:s"=>\$dsh,
				"pid:s"=>\$pid,
				"bid:s"=>\$bid,
				"pop:s"=>\$pop,
				"qtl"=>\$qtl

	-vcf	<file>	input fastq list	
	-anno	<file>	anno
	-out	<dir>	output result dir
	-dsh	<dir>	output work dir
	-pid	<str>	pid
	-bid	<str>	mid
	-pop	<str>	popt
	-qtl

USAGE
	print $usage;
	exit;
}
