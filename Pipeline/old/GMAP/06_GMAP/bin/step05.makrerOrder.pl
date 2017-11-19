#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$out,$dsh,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gen:s"=>\$fIn,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($fIn and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open In,$fIn;
open SH,">$dsh/step05.markerOrder.sh";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($lg,$file)=split(/\s+/,$_);
	if ($popt eq "CP") {
		my $loc="$file";
		my $pwd;
		for (my $i=0;$i<10;$i++) {
			print SH "perl $Bin/bin/calculatePWDviaQsub.pl -i $file -k cycle$i.$lg -d $out/ -n 200 && ";
			print SH "perl $Bin/bin/linkagePhase.pl -i $out/cycle$i.$lg.pwd.detail -g $file -k cycle$i.$lg -d $out/ && ";
			print SH "perl $Bin/bin/extractPwdViaLP.pl -i $out/cycle$i.$lg.pwd.detail -l $out/cycle$i.$lg.loc -k cycle$i.$lg -d $out && ";
			print SH "sgsMap -L $loc -P $pwd -K $out/cycle$i.$lg &&";
			print SH "perl $Bin/bin/smooth-CP.pl -m $out/cycle$i.$lg.map -l $loc -k cycle$i -d $out\n";
		}
	}else{
		for (my $i=0;$i<10;$i++) {
			print SH "MSTMap $file $out/cycle$i.$lg.out &&";
			print SH "perl $Bin/bin/smooth-NOCP.pl -i $file -o $out/cycle$i.$lg.mstmap -t $out/cycle$i.$lg.out -D 0.7 -M 0.7 && ";
			$file="$out/cycle$i.$lg.mstmap.loc";
		}
		print SH "\n";
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
  -gen	<file>	input gen file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -h         Help

USAGE
        print $usage;
        exit;
}
