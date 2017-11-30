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
##############################Mstmap 流程
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($locfile,$dOut,$type,$fKey,$min_theh,$min_miss,$posifile);
GetOptions(
				"help|?" =>\&USAGE,
				"loc:s"=>\$locfile,
				"p:s"=>\$posifile,
				"t:s"=>\$type,
				"k:s"=>\$fKey,
				"o:s"=>\$dOut,
				"min_H:s"=>\$min_theh,
				"min_X:s"=>\$min_miss,
				
				) or &USAGE;
&USAGE unless ($locfile and $type and $fKey and $dOut );

#
# command line options 
#
mkdir ($dOut) if (!-d $dOut);
$locfile = Cwd::abs_path($locfile) ;
$dOut    = Cwd::abs_path($dOut) ;

$min_theh||=0.95;
$min_miss||= 0.9;


my $cur_time;
my $command;
open ( LOG , ">$dOut/MSTmap_pro.log" ) or die $!;

#
# loc format 
#
$type = "RIL".$1 if ($type =~/[Rr][Ii](\d+)/);
$type = "RIL2" if ($type eq "F2") ;
$type = "BC1" if($type=~/BC/i);							#马立祥
#
# MSTmap fails to infer linkage phase for DH, so adjust loc file according to linkage phase for MSTmap
#
#if ($type =~ /[Dd][Hh]/){
if ($type =~ /[Dd][Hh]/ || $type =~/BC1/i){					#马立祥
	$command = "perl $Bin/calculate_pwd_DH.pl -i $locfile -o $dOut/$fKey -t $type";

	$cur_time = &GetTime ;
	print LOG $cur_time,"\n" ;
	print LOG $command,"\n" ;

	`$command `;

	$cur_time = &GetTime ;
	print LOG $cur_time,"\n" ;
}

#if ($type =~ /[Dd][Hh]/){
if ($type =~ /[Dd][Hh]/ || $type =~/BC1/i){					#马立祥
	$locfile = "$dOut/$fKey.DH.loc"; ## important
}

#
# generate MSTmap input 
#
$command = "perl $Bin/loc_to_MSTloc.pl -i $locfile -t $type -k $fKey -o $dOut/$fKey.mst.loc " ;

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";
print LOG "$command\n\n";

`$command `;

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";


#
# marker ordering using MSTmap
#
my $map_dir = "$dOut/MST_ordering";
mkdir $map_dir unless (-d "$map_dir") ;

$command = "$Bin/MSTmap $dOut/$fKey.mst.loc $map_dir/$fKey.mst.map 1>$map_dir/$fKey.mst.o  2>$map_dir/$fKey.mst.e";

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";
print LOG "$command\n\n";

`$command `;

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";


#
# extract MSTmap result 
#
##############liangsh add $posifile and modified the format MST result scripts
if ($posifile) {
	$command = "perl $Bin/format_MST_result_ref.pl -i $map_dir/$fKey.mst.map -l $locfile -m $map_dir/$fKey.mst.o -p $posifile -k $fKey -o $dOut";
}else{
	$command = "perl $Bin/format_MST_result.pl -i $map_dir/$fKey.mst.map -l $locfile -k $fKey -o $dOut";
	}

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";
print LOG "$command\n\n";

`$command `;

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";

#
# evaluate map: heat map and haplotype map  
#
my $draw_heatmap ;
if ( $type =~ /[Dd][Hh]/ || $type =~/BC1/i ) {
	$draw_heatmap = "perl $Bin/draw/draw_heatmap.pl -i $dOut/$fKey.map -p $dOut/$fKey.DH.pwd -k  $fKey -d $dOut/map" ;
}elsif ($posifile) {
	$draw_heatmap = " perl $Bin/draw/draw_MST_heatmap.pl -i $dOut/$fKey.mst_order.o -k  $fKey -d $dOut/map" ;
}else{
	$draw_heatmap = " perl $Bin/draw/draw_MST_heatmap.pl -i $dOut/MST_ordering/$fKey.mst.o -k  $fKey -d $dOut/map" ;
}
my $draw_haplomap = "perl $Bin/draw/draw_haplotype_map.pl  -i $dOut/$fKey.order.loc  -o $dOut/map/$fKey\_haplo  " ;

$cur_time = &GetTime ;
print LOG $cur_time,"\n" ;
print LOG $draw_heatmap,"\n", $draw_haplomap ,"\n";

`$draw_heatmap `;
`$draw_haplomap `;

$cur_time = &GetTime ;
print LOG $cur_time ,"\n";

#
# statistic singleton and suspicious genotypes 
#

$command = "perl $Bin/stat_suspicious_genotype.pl -i $dOut/$fKey.order.loc -o $dOut/$fKey.order.stat";

print LOG "$command\n\n";

`$command `;

#
# imput missing observations and genotyping errors
#

my $smooth_dir = "$dOut/smooth";
mkdir $smooth_dir unless (-d "$smooth_dir");

$command = "perl $Bin/smooth.pl -i $dOut/$fKey.order.loc -o $smooth_dir/$fKey.corr -D $min_theh -M $min_miss";

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";
print LOG "$command\n\n";

`$command `;

$cur_time = &GetTime;
print LOG $cur_time ,"\n\n";

close (LOG) ;
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
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: 单轮 MSTmap+smooth 排图
  Options:
  
  -help|?		USAGE
  -loc			input loc file, forced  
  -k			output file stem, forced
  -o			output directory, forced 
  -t			population type, DH,Ri{2,3,4,5,6...}, forced 
  -min_H		smooth min threshold to change a loc, optional, default [0.95]
  -min_X		smooth min misratio to change a loc, optional, default [0.9]
  
USAGE
	print $usage;
	exit;

}

