#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($lg,$gen,$dOut,$dsh,$ref,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"lg:s"=>\$lg,
	"gen:s"=>\$gen,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dsh,
	"ref:s"=>\$ref,
	"popt:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($lg and $gen and $dOut and $dsh and $popt);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
$lg=ABSOLUTE_DIR($lg);
$gen=ABSOLUTE_DIR($gen);
open SH,">$dsh/step05-1.split.sh";
if ($popt eq "CP") {
	print SH "perl $Bin/bin/splitbyLG-CP.pl -l $lg -i $gen -d $dOut/ -t $popt ";
}else{
	print SH "perl $Bin/bin/splitbyLG-NOCP.pl -l $lg -i $gen -d $dOut/ -t $popt";
}
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/step05-1.split.sh";
print $job;
`$job`;
if ($ref) {
	open SH,">$dsh/step05-2.ref.sh";
	open In,"$dOut/pri.marker.list";
	open Out,">$dOut/ref.marker.list";
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($lg,$file)=split(/\s+/,$_);
		my $map=$file;
		$map=s/\.marker/\.map/g;
		if ($popt eq "CP") {
			print SH "crosslink_group --inp=$file --outbase=$dOut/$lg  --min_lod=-1 --ignore_cxr=1  --randomise_order=1 && ";
			print SH "perl $Bin/bin/smooth-CP.pl -l $dOut/$lg\000.loc -m $dOut/$map -k $lg -d ./\n";
		}else{
			print SH "perl $Bin/bin/smooth-NOCP.pl -i $dOut/$file -o $dOut/$lg.correct.loc -m $dOut/$map";
		}
		print Out "$lg\t$dOut/$lg.correct.loc\n";
	}
	close In;
	close SH;
	close Out;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dsh/step05-2.ref.sh";
	print $job;
	`$job`;
}

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
  -lg	<file>	input lg file name
  -gen	<file>	input gen file
  -out	<dir>	output file dir
  -dsh	<dir>	output worksh dir
  -popt	<srt>	population type
  -ref			modify by ref default off
  -h         Help

USAGE
        print $usage;
        exit;
}
