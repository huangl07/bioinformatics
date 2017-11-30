#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$popt,$out,$bin,$pid,$mid,$nchr,$bin);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"popt:s"=>\$popt,
	"nchr:s"=>\$nchr,
	"out:s"=>\$out,
	"pid:s"=>\$pid,
	"mid:s"=>\$mid,
	"bin:s"=>\$bin,
			) or &USAGE;
&USAGE unless ($vcf and $popt and $out and $pid and $mid);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
my $dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
open LOG,">$out/work_sh/gmap.$BEGIN_TIME.log";
my $step=1;
if ($step == 1) {
	print Log "########################################\n";
	print Log "variant-convert \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.variant-convert.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid  -dsh $dsh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "binner-calculate \n",my $time=time();
	print Log "########################################\n";
	my $marker=ABSOLUTE_DIR("$out/01.vcf--convert/pop.marker");
	my $job="perl $Bin/bin/step02.variant-bin.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid -dsh $dsh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "mlod-calc \n",my $time=time();
	print Log "########################################\n";
	my $marker;
	if ($bin) {
		$marker=ABSOLUTE_DIR("$out/02.vcf--convert/Total.bin.marker");
	}else{
		my $marker=ABSOLUTE_DIR("$out/01.vcf--convert/pop.marker");
	}
	my $job="perl $Bin/bin/step03.mlod-calc.pl -input $marker -out $out/03.mlod-calc -dsh $dsh \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "grouping \n",my $time=time();
	print Log "########################################\n";
	my $mlod=ABSOLUTE_DIR("$out/03.mlod-calc/Total.mlod");
	my $marker;
	if ($bin) {
		$marker=ABSOLUTE_DIR("$out/02.vcf--convert/Total.bin.marker");
	}else{
		my $marker=ABSOLUTE_DIR("$out/01.vcf--convert/pop.marker");
	}
	my $job="perl $Bin/bin/step04.grouping.pl -marker $marker -popt $popt -out $out/04.grouping -nchr $nchr -dsh $dsh ";
	if ($ref ) {
		$jobs.= "-ref";
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 5) {
	print Log "########################################\n";
	print Log "map cycle1 \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.markerOrder.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 6) {
	print Log "########################################\n";
	print Log "map cycle2 \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.markerOrder.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 7) {
	print Log "########################################\n";
	print Log "map cycle3 \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.markerOrder.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 8) {
	print Log "########################################\n";
	print Log "map cycle8 \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.markerOrder.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 9) {
	print Log "########################################\n";
	print Log "map cycle9 \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step05.markerOrder.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}

if ($step == 10) {
	print Log "########################################\n";
	print Log "evalutaion\n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step06.mapEvaluation.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid \n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}

close LOG;
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
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -vcf	<file>	input file name
  -out	<dir>	output dir of filename
  -popt	<str>	population type CP/BCi/Fi/Rix/
  -bin		binmap or not 
  -h         Help

USAGE
        print $usage;
        exit;
}
