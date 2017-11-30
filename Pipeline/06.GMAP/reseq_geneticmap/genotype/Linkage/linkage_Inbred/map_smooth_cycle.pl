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
my ($locfile,$type,$fKey,$dOut,$record,$MSTmap,$min_theh,$min_miss,$nCycle,$posifile);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$locfile,
				"o:s"=>\$dOut,
				"k:s"=>\$fKey,
				"t:s"=>\$type,
				"p:s"=>\$posifile,
				
				"record"=>\$record,
				"MSTmap"=>\$MSTmap,
		
				"n:s"=>\$nCycle,
				## smooth
				"min_H:s"=>\$min_theh,
				"min_X:s"=>\$min_miss,
				
				) or &USAGE;
&USAGE unless ($locfile and $type and $fKey and $dOut);
&USAGE unless ( $record || $MSTmap );

mkdir ($dOut) if (!-d $dOut);
$dOut = Cwd::abs_path($dOut);
$locfile = Cwd::abs_path($locfile);

#
# command line option 
#
$nCycle ||= 10 ;
$min_theh||=0.95;
$min_miss||= 0.8;


my $cur_time ;
my $command ;
open (Log , ">$dOut/$fKey.mapcycle.log") or die $!;
my $key = $fKey;

#
# MSTmap
#
if ($MSTmap){

	$fKey = $key."\_MST";
	my $mstOut = "$dOut/MSTmap/";
	mkdir ($mstOut) if (!-d $mstOut);
	
	#
	# first cycle 
	#
	my $cycle = 0;
	###########################liansh add $posifile

	$command = "perl $Bin/MSTmap_smooth_main.pl -loc $locfile -k $fKey\_cycle$cycle -t $type -o $mstOut/cycle_$cycle -min_H $min_theh -min_X $min_miss" ;
	if ($posifile) {
		$command .= " -p $posifile ";
	}

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
	print Log $command,"\n" ;

	`$command `;

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;


	#
	# MST - smooth cycle 
	#
	while (1) {
		my %misratio ;
		my %singlratio;
	
		#
		# statistic and test convergence 
		#
		open (Xls, "$mstOut/cycle_$cycle/smooth/$fKey\_cycle$cycle.corr_Diff.xls" ) or die $!;
		while (<Xls>) {
			s/\%//g;
			next if (/^\s*$/) ;
			my @temp =split ;
			if (/before/) {
				$misratio{before} = $temp[4] ;
				$singlratio{before} = $temp[5] ;
			}elsif (/after/){
				$misratio{after} = $temp[4] ;
				$singlratio{after} = $temp[5] ;
			}
		}
		last if (($misratio{before} < 3 || $misratio{before} - $misratio{after} < 0.02 ) 
				&& $singlratio{before} - $singlratio{after} < 0.05 ) ;
		close (Xls) ;
		

		my $pre_cycle = $cycle;
		$cycle++ ;
		last if( $cycle > $nCycle);

		#
		# next cycle 
		#
		$command = "perl $Bin/MSTmap_smooth_main.pl -loc $mstOut/cycle_$pre_cycle/smooth/$fKey\_cycle$pre_cycle.corr.loc -k $fKey\_cycle$cycle -t $type -o $mstOut/cycle_$cycle  -min_H $min_theh -min_X $min_miss";

		if ($posifile) {
			$command .= " -p $posifile";
		}

		$cur_time = &GetTime ;
		print Log $cur_time,"\n" ;
		print Log $command, "\n" ;

		`$command ` ;

		$cur_time = &GetTime ;
		print Log $cur_time,"\n" ;
		
	}


	#
	# statistic MST-smooth cycle info 
	#
	$command =  "perl $Bin/stat_cycle_info.pl  -i $mstOut -k $fKey  -o $mstOut";

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
	print Log $command,"\n" ;

	`$command ` ;

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;

}

##############################################################################################################################
#
# RECORD
#
if ($record) {
	
	$fKey = $key."\_RECORD";
	my $recordOut = "$dOut/record";
	mkdir ($recordOut) if (!-d $recordOut);


	#
	# first cycle 
	#

	my $cycle = 0;
	$command = "perl $Bin/RECORD_smooth_main.pl -loc $locfile -k $fKey\_cycle$cycle -t $type -o $recordOut/cycle_$cycle -min_H $min_theh -min_X $min_miss " ;
	$command .= " -p $posifile" if ($posifile);

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
	print Log $command,"\n" ;

	`$command `;

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;


	#
	# RECORD - smooth cycle 
	#
	while (1) {
		my %misratio;
		my %singlratio;
		
		#
		# stat and test convergence
		#
		open (Xls, "$recordOut/cycle_$cycle/smooth/$fKey\_cycle$cycle.corr_Diff.xls" ) or die $!;
		while (<Xls>) {
			s/\%//g;
			next if (/^\s*$/) ;
			my @temp =split ;
			if (/before/) {
				$misratio{before} = $temp[4] ;
				$singlratio{before} = $temp[5] ;
			}elsif (/after/){
				$misratio{after} = $temp[4] ;
				$singlratio{after} = $temp[5] ;
			}
		}

		last if (($misratio{before} < 2 || $misratio{before} - $misratio{after} < 0.1) 
				&& $singlratio{before} - $singlratio{after} < 0.05 ) ;
		close (Xls) ;

		my $pre_cycle = $cycle;
		$cycle++ ;
		last if( $cycle > $nCycle);
		
		#
		# next cycle 
		#
		$command = "perl $Bin/RECORD_smooth_main.pl -loc $recordOut/cycle_$pre_cycle/smooth/$fKey\_cycle$pre_cycle.corr.loc -k $fKey\_cycle$cycle -t $type -o $recordOut/cycle_$cycle  -min_H $min_theh -min_X $min_miss " ;
		$command .= " -p $posifile" if ($posifile);
		$cur_time = &GetTime ;
		print Log $cur_time,"\n" ;
		print Log $command,"\n" ;

		`$command ` ;

		$cur_time = &GetTime ;
		print Log $cur_time,"\n" ;
		
	}

	#
	# statistic RECORD - smooth cycle 
	#
	
	$command =  "perl $Bin/stat_cycle_info.pl  -i $recordOut -k  $fKey  -o $recordOut";

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
	print Log $command,"\n" ;

	`$command ` ;

	$cur_time = &GetTime ;
	print Log $cur_time,"\n" ;
}



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

Usage: 近交群体多轮纠错排图流程， 排图使用MSTmap 或 RECORD排序算法
  Options:

  -i			input loc file, forced
  -o			output directory, forced 
  -k			output file stem, forced  
  -t			population type, forced 
  
  -record		record mapping method
  -MSTmap		MST mapping method 

  -n			maximum cycle numbers allowed, default [10] 
  -min_H		smooth min threshold to change a loc 0.95
  -min_X		smooth min ratio to change a loc 0.8
 
  -help			USAGE

USAGE
	print $usage;
	exit;
}

