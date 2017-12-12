#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$popt,$out,$bin,$pid,$mid,$nchr,$ref,$step);
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
	"bin"=>\$bin,
	"ref"=>\$ref,
	"step:s"=>\$step,
			) or &USAGE;
&USAGE unless ($vcf and $popt and $out and $pid and $mid and $nchr);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
my $dsh="$out/work_sh";
mkdir $dsh if (!-d $dsh);
open Log,">$out/work_sh/gmap.$BEGIN_TIME.log";
$step||=1;
if ($step == 1) {
	print Log "########################################\n";
	print Log "variant-convert \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.variant-convert.pl -vcf $vcf -popt $popt -out $out/01.vcf-convert -pid $pid -mid $mid -dsh $dsh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
	$step++ if (!$bin);
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "binner-calculate \n",my $time=time();
	print Log "########################################\n";
	my $marker=ABSOLUTE_DIR("$out/01.vcf-convert/pop.filtered.marker");
	my $job="perl $Bin/bin/step02.variant-bin.pl -marker $marker -popt $popt -out $out/02.binner -dsh $dsh\n";
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
		$marker=ABSOLUTE_DIR("$out/02.binner/Total.bin.marker");
	}else{
		$marker=ABSOLUTE_DIR("$out/01.vcf-convert/pop.filtered.marker");
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
		$marker=ABSOLUTE_DIR("$out/02.binner/Total.bin.marker");
	}else{
		$marker=ABSOLUTE_DIR("$out/01.vcf-convert/pop.filtered.marker");
	}
	my $job="perl $Bin/bin/step04.grouping.pl -marker $marker -mlod $mlod -out $out/04.grouping -nchr $nchr -popt $popt -dsh $dsh ";
	if ($ref ) {
		$job.= "-ref";
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
	my $lg=ABSOLUTE_DIR("$out/04.grouping/Total.lg");
	my $marker;
	if ($bin) {
		$marker=ABSOLUTE_DIR("$out/02.binner/Total.bin.marker");
	}else{
		$marker=ABSOLUTE_DIR("$out/01.vcf-convert/pop.filtered.marker");
	}
	my $job="perl $Bin/bin/step05.markerOrder.pl -lg $lg -gen $marker -popt $popt -out $out/05.map-cycle1 -dsh $dsh  -cycle 1\n";
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
	my $gen=ABSOLUTE_DIR("$out/05.map-cycle1/marker.list");
	my $job="perl $Bin/bin/step05.markerOrder.pl -gen $gen -popt $popt -out $out/06.map-cycle2 -dsh $dsh -cycle 2\n";
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
	my $gen=ABSOLUTE_DIR("$out/06.map-cycle2/marker.list");
	my $job="perl $Bin/bin/step05.markerOrder.pl -gen $gen -popt $popt -out $out/07.map-cycle3 -dsh $dsh -cycle 3\n";
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
	print Log "map cycle4 \n",my $time=time();
	print Log "########################################\n";
	my $gen=ABSOLUTE_DIR("$out/07.map-cycle3/marker.list");
	my $job="perl $Bin/bin/step05.markerOrder.pl -gen $gen -popt $popt -out $out/08.map-cycle4 -dsh $dsh -cycle 4\n";
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
	print Log "map cycle5 \n",my $time=time();
	print Log "########################################\n";
	my $gen=ABSOLUTE_DIR("$out/08.map-cycle4/marker.list");
	my $job="perl $Bin/bin/step05.markerOrder.pl -gen $gen -popt $popt -out $out/09.map-cycle5 -dsh $dsh -cycle 5\n";
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
	my $marker=ABSOLUTE_DIR("$out/01.vcf-convert/pop.filtered.marker");
	my $dmap=ABSOLUTE_DIR("$out/09.map-cycle5");
	my $job="perl $Bin/bin/step06.mapEvaluation.pl -dmap $dmap -popt $popt -out $out/10.mapEvalue --mark $marker -dsh $dsh ";
	if ($ref) {
		$job.="-ref \n";
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}

close Log;
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
  -nchr	<num>	number of chromosome number
  -pid	<str>	paternal id
  -mid	<str>	maternal id
  -bin		binmap or not 
  -ref		ref base or not
  -h         Help

USAGE
        print $usage;
        exit;
}
