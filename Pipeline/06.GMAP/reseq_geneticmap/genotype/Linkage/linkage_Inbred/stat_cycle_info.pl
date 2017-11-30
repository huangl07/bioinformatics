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
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dIn,$dOut,$fKey);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$dIn,
				"k:s"=>\$fKey,
				"o:s"=>\$dOut,
				) or &USAGE;
&USAGE unless ($dIn and $dOut);

mkdir ( $dOut ) if (! -d $dOut) ;
$dIn  = Cwd::abs_path($dIn) ;
$dOut = Cwd::abs_path($dOut) ;

#
# input 
#
my %info_file ;
opendir ( Dr , $dIn) or die $!;
while (defined(my $file = readdir(Dr))) {
	next if ($file !~ /^cycle_(\d+)$/) ;

	$info_file {$1} {"loc"} = "$dIn/cycle_$1/smooth/$fKey\_cycle$1.corr_Diff.xls"  if (-f "$dIn/cycle_$1/smooth/$fKey\_cycle$1.corr_Diff.xls") ;
	$info_file {$1} {"map"} = "$dIn/cycle_$1/$fKey\_cycle$1.map"  if (-f "$dIn/cycle_$1/$fKey\_cycle$1.map") ;
	$info_file {$1} {"noise"} = "$dIn/cycle_$1/$fKey\_cycle$1.order.stat.xls"  if (-f "$dIn/cycle_$1/$fKey\_cycle$1.order.stat.xls") ;
}
closedir(Dr) ;

#
# statistic and output 
#
open (STAT , ">", "$dOut/$fKey\_cycle_info.xls") or die $!;
print STAT "cycle\tnloc\tnind\tsum\tmiss\tsingleton\tnoise\tmiss\%\tsingleton\%\tnoise\%\tgap_5\tgap_5\%\tmapdistance\t\n" ;

my ($loc,$ind,$sum,$miss,$singleton,$noise,$miss_ratio,$singleton_ratio,$noise_ratio,$map_dist);
foreach my $cycle_num ( sort {$a <=> $b} keys %info_file ) {

	#
	# miss and singleton 
	#
	open (LOC , $info_file{$cycle_num}{"loc"} ) or die $!;
	while (<LOC>) {
		next if (!/before/) ;
		(undef,$sum,$miss,$singleton,$miss_ratio,$singleton_ratio) = split ;
	}
	close (LOC) ;
	
	#
	# noise
	#
	my @gap = ();
	my $distance = 0 ;
	open (Noise , $info_file{$cycle_num}{"noise"} ) or die $!;
	while (<Noise>) {
		next if (/^\s*$/ || /^\#/) ;

		if (/nloc\D*(\d+)$/) {
			$loc = $1 ;
			next ;
		}elsif (/nind\D*(\d+)/){
			$ind = $1 ;
			next ;
		}
		($sum , $miss , $noise , $miss_ratio ,$noise_ratio ) = split ;
	}
	close (Noise) ;
	
	#
	# map info 
	#
	my $map_dist = 0 ;
	my $i = 0 ;
	open (Map , $info_file{$cycle_num}{"map"} ) or die $!;
	while (<Map>) {

		next if (/^group/ || /^;/ || /^$/) ;
		
		my ($marker, $cm) = split ;
		$map_dist = $cm if ( $cm > $map_dist ) ;
		push @gap ,$cm ;
	}
	close (Map) ;
	
	my $gap_5 = 0 ;
	map {$gap_5 ++ if ( ($gap[$_] - $gap[$_-1]) < 5) } 1..@gap - 1 ;

	my $gap_5_ra = $gap_5 /($loc -1 )*100 ;
	$gap_5_ra .= "%" ;

	print STAT join( "\t", 
		"cycle$cycle_num",
		$loc,
		$ind,
		$sum,
		$miss,
		$singleton,
		$noise,
		$miss_ratio,
		$singleton_ratio,
		$noise_ratio, 
		$gap_5, 
		$gap_5_ra,
		$map_dist
	),"\n";
}

close (STAT) ;
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

Usage: 统计多轮纠错结果
  Options:
  -i			input dir(MSTmap and RECORD analysis), forced 
  -k			output file stem, forced 
  -o			output directory 
  -h			usage

USAGE
	print $usage;
	exit;
}
