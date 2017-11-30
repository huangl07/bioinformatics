#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$dOut,$log,$onlyFirst,$start,$end,$stepSize,$minGroup,$maxGroup,$nChro,$minMarkerNum,@step);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				"n:s"=>\$nChro,

				"p:i"=>\$onlyFirst,
				"m:i"=>\$minMarkerNum,

				"b:s"=>\$start,
				"e:s"=>\$end,
				"s:s"=>\$stepSize,
				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,
				"step{,}:s"=>\@step

				) or &USAGE;
&USAGE unless ($fIn and $fKey and $nChro);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$dOut = Cwd::abs_path($dOut);
$fIn = Cwd::abs_path($fIn);

$minMarkerNum||=400;
$onlyFirst||=200;

$start||=3;        ## lod start 
$end||=20;         ## lod end 
$stepSize||=1;     ## lod step size
$minGroup||=20;    ## minimum threshold for group size
$maxGroup||=1000;   ## maximum threshold for group size
@step||=qw(10,100);
#-------------------------------------------------------------------
# Process
#-------------------------------------------------------------------

### create makefile 
open SH,">$dOut/work.sh";
while ($minGroup < $maxGroup) {
	print SH "perl $Bin/determine_linkage_group.pl -i $dOut/$fKey.mLOD -k $fKey -d $dOut -n $nChro -b $start -e $end -s $stepSize -minGroup $minGroup -maxGroup $maxGroup &&\n";
	$minGroup+=$step[0];
	$maxGroup+=$step[1];
}
close SH;
open (MAKE,">$dOut/Makefile") or die $!;

print MAKE "#------------ determine linkage groups -------------------\n";
print MAKE "$dOut/$fKey.*.lg : $dOut/$fKey.mLOD \\\n $dOut/calculate_mLOD.check \n";
print MAKE "\tsh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $dOut/work.sh \n";

print MAKE "#------------ calculate modified LOD score ----------------\n";
print MAKE "$dOut/$fKey.mLOD : $fIn\n";
print MAKE "\tperl $Bin/calculate_mLOD_via_qsub.pl -i $fIn -d $dOut -k $fKey -p $onlyFirst -m $minMarkerNum 1>>$dOut/make.o 2>>$dOut/make.e \n";


print MAKE "#-----------------------\n";
print MAKE "clean :\n";
print MAKE "\t\@echo \"collapse project\"\n";
print MAKE "\t-rm -r $dOut/*.log $dOut/*.check $dOut/*.mLOD $dOut/*.hash $dOut/all_scheme $dOut/calculate_mLOD_dir $dOut/make.o $dOut/make.e \n";
print MAKE "\t\@echo \"collapse project done\"\n";
print MAKE ".PHONY : clean";

close (MAKE) ;

### make 

`make -C $dOut`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:
		A small pipeline for linkage grouping.
Usage:
  Options:
  -i	<file>	Genotype file, BMK format, forced
  -k	<str>	Key of output file,forced
  -d	<str>	Directory where output file produced,optional,default [./]
  -n	<int>   Species\' chromosome number,forced
  
  -p	<int>	Integer skip number marker at split file, default 200
  -m	<int>	Minimum threshold of markers\' number for qsub, optional, default [400]

  -b	<int>   Start lod,optional,default [3]
  -e	<int>   End lod,optional,default [20]
  -s	<int>   Stepsize of lod,optional,default [1]

  -minGroup     <int>   Minimum threshold of group size,optional,default [20]
  -maxGroup     <int>   Maximum threshold of group size,optional,default [1000]
  -step	<int>	step to adjust minGroup and maxGroup 
				default(10,100) means determinGroup by several parameters:
				20	1000
				30	900
				40	800
				50	700
				60	600
				70	500
				80	400
				90	300
				100	200
				
  -h		Help

USAGE
	print $usage;
	exit;
}
