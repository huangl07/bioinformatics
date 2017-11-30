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
# ==============================================================
# Get Options
# ==============================================================
my ($fIn,$fAlign,$popt,$nChro,$fKey,$dOut,$fScore,$log, $modify, $minGroup, $maxGroup, $step, $only);

GetOptions(
		"help|?" =>\&USAGE,
		"i:s"=>\$fIn,	
		"k:s"=>\$fKey,
		"d:s"=>\$dOut,
		
		"popt:s"=>\$popt,
		"nchro:i"=>\$nChro,
		"a:s"=>\$fAlign,
		"modify:i"=>\$modify,
		"s:s"=>\$fScore,

		"minGroup:i"=>\$minGroup,
		"maxGroup:i"=>\$maxGroup,

		"step:s"=>\$step,
		"only"=>\$only,
		) or &USAGE;
&USAGE unless ($fIn and $fKey  and $nChro and $popt);


#===============================================================
# Default optional value 
#===============================================================
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$dOut = Cwd::abs_path($dOut);
$fAlign=Cwd::abs_path($fAlign) if (defined $fAlign and -f $fAlign);
$minGroup ||= 20;
$maxGroup ||= 500;
$step ||=1;
$modify = 1 unless (defined $modify);

#
# Output log file 
#
open ($log,">","$dOut/$fKey.linkage.".time.".log") or die $!;
my $startTime =  GetTime();
my $user      =  `whoami`;chomp $user;
my $workDir   =  `pwd`;chomp $workDir;
my $task      =  "perl $Bin/$Script ".join(" ",@Original_ARGV);

print $log "######################################\n";
print $log "$startTime\n";
print $log "user:	$user\n";
print $log "work directory:	$workDir\n";
print $log "job:	$task\n";
print $log "######################################\n\n";

print $log "input genotype file:	$fIn\n";
print $log "analysis directory:	$dOut\n";
print $log "project name:	$fKey\n";
print $log "population type:	$popt\n";
print $log "nr. of chromosome:	$nChro\n";

print $log "minimum group size: $minGroup\n";
print $log "maximum group size: $maxGroup\n";

print $log "step of pipeline:	$step\n";
print $log "\n\n";

#===============================================================
# pipeline
#===============================================================



	print "--- perform linkage analysis ---\n";
	print $log "--- perform linkage analysis ---\n";

	my $job = "";
	if ($popt eq "CP") {
		$job  = "perl $Bin/bin/linkage_CP/linkage_CP.pl -i $fIn -k $fKey -d $dOut -n $nChro  -minGroup $minGroup -maxGroup $maxGroup -step $step -modify $modify" ;
		$job .= " -ref $fAlign " if (defined $fAlign and -f $fAlign);
		$job .= " -only" if ($only) ;
	}#elsif ($popt eq "F2" || $popt =~ /Ri/ || $popt eq "BC1" || $popt eq "DH") {
	elsif ($popt eq "F2" || $popt =~ /Ri/ || $popt =~/BC\d+/ || $popt eq "DH"){
#		$job  = "perl $Bin/bin/linkage_Inbred/linkage_Inbred.pl -i $fIn -d $dOut -k $fKey -popt $popt -nChro $nChro -record -MSTmap -minGroup $minGroup -maxGroup $maxGroup -step $step -modify $modify";
		$job  = "perl $Bin/bin/linkage_Inbred/linkage_Inbred.pl -i $fIn -d $dOut -k $fKey -popt $popt -nChro $nChro -MSTmap -minGroup $minGroup -maxGroup $maxGroup -modify $modify";
		$job .= " -ref $fAlign " if (defined $fAlign and -f $fAlign);
		$job .= " -only" if ($only) ;
	}else{
		die "unfortunately,we can't not deal this $popt Population!\nDo use Joinmap4.1 and Good Luck!\n";
	}

	print "$job\n";
	print $log "$job\n";
	`$job `;
	
	if ($fScore) {
		$job="perl $Bin/bin/stat_mappingSLAF.pl -l $dOut/grouping/*.lg -s $fScore/*.score -o $dOut/Linkage/QC";
		print "$job\n";
		print $log "$job\n";
		`$job `;
	}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close ($log) ;
# ==============================================================
# sub function
# ==============================================================


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

	Usage:
	  Options:
	  -i		<file>	Input genotype file, forced
	  -k		<str>	Key of output file,forced
	  -d		<dir>	Directory where output file produced,optional,default [./]
	  -popt		<str>	Population type, forced
	  -s		<file>	Score file, optional 
	  -nchro	<int>	The number of chromosome, forced
	  -a		<file>	Alignment file used in grouping with ref, table formated, [MarkerID Chr Start END]
		-modify	<int>   Modify linkage grouping and genotype, 1 or 0, default [1]
	  -minGroup	<int>	Minimum group size, optional, default [20]
	  -maxGroup	<int>	Maximum group size, optional, default [500]
	  
	  -step		<int>	pipiline step 
	  -only             Do only one step 

	  -h			Help

USAGE
	print $usage;
	exit;
}

