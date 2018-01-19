#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$out,$dsh,$popt,$cycle,$lg,$ref);
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
	"ref"=>\$ref,
	"cycle:s"=>\$cycle,
			) or &USAGE;
&USAGE unless ($fIn and $out and $dsh and $popt);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$fIn=ABSOLUTE_DIR($fIn);
mkdir "$out/pri-pwd" if (!-d "$out/pri-pwd");
mkdir "$out/pwd" if (!-d "$out/pwd");

$cycle||=1;
open SH,">$dsh/step05-$cycle.markerOrder.sh";
open List,">$out/marker.list";
open In,$fIn;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($lg,$file)=split(/\s+/,$_);
	if ($popt eq "CP") {
		print SH "crosslink_group --inp=$file --outbase=$out/$lg --mapbase=$out/$lg --min_lod=-1 --randomise_order=1 && ";
		print SH "perl $Bin/bin/smooth-CP.pl -m $out/$lg.sexAver.map -l $out/$lg.loc -k $lg -d $out\n";
	}else{
		print SH "MSTmap $file $out/$lg.out &&";
		print SH "perl $Bin/bin/smooth-NOCP.pl -i $file -m $out/$lg.out -o $out/$lg.correct.marker \n ";
	}
	print List "$lg\t$out/$lg.correct.marker\n";
}
close SH;
close In;
close List;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=10"."G --CPU 1 $dsh/step05-$cycle.markerOrder.sh";
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
  -gen	<file>	input marker list file
  -out	<dir>	output dir
  -dsh	<dir>	worksh dir
  -popt	<srt>	population type
  -h         Help

USAGE
        print $usage;
        exit;
}
