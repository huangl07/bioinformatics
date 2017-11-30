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
#################################loc 
# ------------------------------------------------------------------
# Get Options
# ------------------------------------------------------------------
my ($locfile,$type,$fKey,$fOut,$min_theh,$min_miss,$posifile);
GetOptions(
				"help|?" =>\&USAGE,
				"loc:s"=>\$locfile,
				"o:s"=>\$fOut,
				"k:s"=>\$fKey,
				"t:s"=>\$type,
				"p:s"=>\$posifile,

				"min_H:s"=>\$min_theh,
				"min_X:s"=>\$min_miss,
				) or &USAGE;
&USAGE unless ($locfile and $type and $fKey and $fOut);
#
# command line options 
#
mkdir ($fOut) if (!-d $fOut) ;
$locfile = Cwd::abs_path($locfile);
$fOut    = Cwd::abs_path($fOut);

$min_theh||=0.90;
$min_miss||= 0.8;

my $cur_time ;
open (Log , ">$fOut/$fKey\_pro.log") or die $!;

#
# calculate pair wise data 
#
my $command ;
if ($type !~/[Dd][Hh]/){
	$command = "perl $Bin/calculate_pwd_by_qsub.pl -i $locfile  -d $fOut/ -type $type -k $fKey  -p 100 -m 300" ;
}else {
	$command = "perl $Bin/calculate_pwd_DH.pl -i $locfile  -o $fOut/$fKey";
}


$cur_time = &GetTime ;
print Log $cur_time,"\n" ;
print Log $command,"\n" ;

`$command `;

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;


#
# marker Ordering using RECORD
# 

my $pwdfile = "$fOut/$fKey.pwd" ;
if ($type =~ /[Dd][Hh]/){

	$locfile = "$fOut/$fKey.DH.loc" ;
	$pwdfile = "$fOut/$fKey.DH.pwd" ;

}
$command = "perl $Bin/RECORD.pl -loc $locfile  -pwd $pwdfile -o $fOut/$fKey" ;

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;
print Log $command,"\n" ;

`$command `;

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;

#####order mapfile and locfiel by ref when #######added by liangsh
if ($posifile) {
	$command = "perl $Bin/Record_order_byref.pl -m $fOut/$fKey.map -l $fOut/$fKey.order.loc -p $posifile -o $fOut/$fKey";
	$cur_time = &GetTime ;

	print Log $cur_time,"\n" ;
	print Log $command,"\n" ;
	`$command`;
	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
}

#
# draw heat map and haplotype map
#

my $draw_heatmap = "perl $Bin/draw/draw_heatmap.pl -i $fOut/$fKey.map -p $pwdfile -k  $fKey -d $fOut/map";
my $draw_haplomap = "perl $Bin/draw/draw_haplotype_map.pl -i $fOut/$fKey.order.loc -o $fOut/map/$fKey\_haplo";

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;
print Log $draw_heatmap,"\n", $draw_haplomap ,"\n";

`$draw_heatmap `;
`$draw_haplomap `;

$cur_time = &GetTime ;
print Log $cur_time ,"\n";

#
# stat suspicious genotypes, singletons etc
#

$command = "perl $Bin/stat_suspicious_genotype.pl -i $fOut/$fKey.order.loc -o $fOut/$fKey.order.stat";

print Log "$command\n\n";

`$command `;

#
# smooth 
#
my $smooth_dir = "$fOut/smooth";
mkdir $smooth_dir unless (-d "$smooth_dir");

$command = "perl $Bin/smooth.pl  -i $fOut/$fKey.order.loc -o $smooth_dir/$fKey.corr -D $min_theh -M $min_miss " ;

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;
print Log $command,"\n" ;

`$command `;

$cur_time = &GetTime ;
print Log $cur_time,"\n" ;

close (Log) ;
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
Program:
Version: $version
Contact: Wangml <wangml\@biomarker.com.cn> 

Usage: µ•¬÷ RECORD+smooth ≈≈Õº 
  Options:
  -help			USAGE
  -loc			input loc file, forced 
  -o			output directory, forced 
  -k			output file stem, forced 
  -t			population type, {BC1,Ri+,F2}, forced 
  -p			posifile name, not forced
  -min_H		smooth min threshold to change a loc, optional, default [0.9]
  -min_X		smooth min misratio to change a loc, optional, default [0.8]
 
USAGE
	print $usage;
	exit;
}

