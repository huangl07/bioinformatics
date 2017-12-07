#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$out,$dsh,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"map:s"=>\$dmap,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($dmap and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open In,$fIn;
open SH,">$dsh/step06.mapEvaluation.sh";
if ($popt ne "CP") {
	print SH "perl $Bin/bin/MapMergeNOCP.pl -dmap $dmap -o $out && ";
	print SH "perl $Bin/bin/map2rqtl.pl -l $out/total.marker -m $out/total.map -o $out/total.csvr";
	print SH "Rscript $Bin/bin/emap-NOCP.R --mark $out/total.csvr --out $out --pop $popt\n";
}else{

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
  -gen	<file>	input gen file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -h         Help

USAGE
        print $usage;
        exit;
}
