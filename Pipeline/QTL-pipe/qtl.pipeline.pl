#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($out,$ann,$pop,$btl,$vcf,$trit,$dmap,$step);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"trit:s"=>\$trit,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"vcf:s"=>\$vcf,
	"pop:s"=>\$pop,
	"btl:s"=>\$btl,
	"step:s"=>\$step,
	) or &USAGE;
&USAGE unless ($out and $ann and $trit and $dmap and $vcf and $pop);

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        fq thanslate to fa format
        eg:
        perl $Script -i -o -k -c
Usage:
  Options:
  -dmap	<dir>	input the gmap's result dir(*.csv *.map)
  -trit	<file>	input the qtl's trit list file
  -out	<dir>	output dir
  -ann	<file>	input ann file
  -vcf	<file>  input vcf file
  -pop	<str>	pop type
  -btl		binary trt or not
  -step	<num>	1,2,3
  -h         Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
$pop||="F2";
$step||=1;
$dmap=ABSOLUTE_DIR($dmap);
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$ann=ABSOLUTE_DIR($ann);
$trit=ABSOLUTE_DIR($trit);
mkdir "$out/worksh" if (!-d "$out/worksh");
open Log,">$out/worksh/qtl.$BEGIN_TIME.log";
my $dsh="$out/worksh";
if ($step==1){
	print Log "########################################\n";
	print Log "qtl  \n",my $time=time();
	print Log "########################################\n";
	mkdir "$out/01.qtl" if (!-d "$out/01.qtl");
	my $job="perl $Bin/bin/step01.qtl.pl -dmap $dmap -trit $trit -popt $pop -out $out/01.qtl -dsh $dsh ";
	if($btl){
		$job.="$btl  ";
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	$step++ ;
}
if ($step==2){
	print Log "########################################\n";
	print Log "anno \n",my $time=time();
	mkdir "$out/02.anno" if (!-d "$out/02.anno");
	my $job="perl $Bin/bin/step02.anno.pl -vcf $vcf -ann $ann -qtl $out/01.qtl -out $out/02.anno -dsh $dsh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	$step++ ;
}
if ($step==3){
	print Log "########################################\n";
	print Log "merge \n",my $time=time();
	mkdir "$out/03.report" if (!-d "$out/03.report");
	my $job="perl $Bin/bin/step03.merge.pl -in $out -out $out/03.report -pop $pop -dsh $dsh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
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
