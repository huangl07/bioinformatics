#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$pid,$bid,$popt,$ann);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"pid:s"=>\$pid,
	"bid:s"=>\$bid,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($vcf and $out and $bid and $ann);
mkdir $out if (!-d $out);
$popt||="F2";
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$ann=ABSOLUTE_DIR($ann);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
open SH,">$out/work_sh/bsa.sh";
my @bid=split(/\,/,$bid);
if (scalar @bid == 1) {
	print SH "perl $Bin/bin/bin/mutmap.pl -vcf $vcf -out $out/index-calc.result -bid $bid ";
	if ($pid) {
		print SH "-pid $pid && ";
	}else{
		print SH "&& ";
	}
	print SH "Rscript $Bin/bin/bin/slidingwin.R --infile $out/index-calc.result --outfile $out/sliding-win --col 1,2,10 --win 2000000 --step 10000 --method bp && ";
	print SH "Rscript $Bin/bin/bin/manhattan.R --input $out/sliding-win.result --output $out/bsa --col 1,3,4,5 && ";
	print SH "perl $Bin/bin/bin/region-variant.pl -i $out/index-calc.result -o $out/region.threshold.out -r $out/sliding-win.threshold.select\n";
	print SH "perl $Bin/bin/bin/region-variant.pl -i $out/index-calc.result -o $out/region.fdr.out -r $out/sliding-win.fdr.select\n";
}else{
	print SH "perl $Bin/bin/qtlseq.pl -vcf $vcf -out $out/index-calc.result -bid $bid -popt $popt ";
	if ($pid) {
		print SH "-pid $pid && ";
	}else{
		print SH "&& ";
	}
	print SH "Rscript $Bin/bin/slidingwin.R --infile $out/index-calc.result --outfile $out/sliding-win --col 1,2,14,15,16 --win 2000000 --step 10000 method bp && ";
	print SH "Rscript $Bin/bin/manhattan.R --infile $out/sliding-win.result --outfile $out/bsa --col 1,3,4,5,6 && ";
	print SH "perl $Bin/bin/region-variant.pl -i $out/sliding-win.threshold.select -o $out/region.out -a $out/index-calc.result";
}

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
  -vcf	<file>	input vcf file
  -out	<dir>	output dir
  -pid	<id>	input p id
  -bid	<id>	input b id
  -ann	<file>	input ann file
  -h         Help

USAGE
        print $usage;
        exit;
}
