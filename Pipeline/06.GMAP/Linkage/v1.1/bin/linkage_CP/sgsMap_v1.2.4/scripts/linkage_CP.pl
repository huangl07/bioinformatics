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
my ($fIn,$fKey,$dOut,$minGroup,$maxGroup,$nchro,$cycle_t,$step,$only,$log);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				
				"n:s"=>\$nchro,
				
				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,
				
				"cycle:s"=>\$cycle_t,

				"step:s"=>\$step,
				"only"=>\$only,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $nchro);
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$dOut = Cwd::abs_path($dOut);

$step||= 1;
$minGroup||=20;
$maxGroup||=500;
$cycle_t||=5;

open ($log,">","$fKey.log") or die $!;
&printLogHead($log);

#-------------------------------------------------------------------
# Pipeline
#-------------------------------------------------------------------

# step 1 : determine linkage group 
if ($step == 1) {
	print "determine linkage groups\n";
	print $log "determine linkage groups\n";
	
	my $groupDir = "$dOut/grouping";

	print "perl $Bin/grouping/grouping.pipeline.pl -i $fIn -k $fKey -d $groupDir -n $nchro -minGroup $minGroup -maxGroup $maxGroup\n";
	print $log "perl $Bin/grouping/grouping.pipeline.pl -i $fIn -k $fKey -d $groupDir -n $nchro -minGroup $minGroup -maxGroup $maxGroup\n";
	`perl $Bin/grouping/grouping.pipeline.pl -i $fIn -k $fKey -d $groupDir -n $nchro -minGroup $minGroup -maxGroup $maxGroup`;
	
	die "determine linkage group failed" unless (glob("$dOut/grouping/*.lg")) ;

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
}

# step 2 : construct linkage map each linkage group 
if ($step == 2) {
	print "mapping process\n";
	print $log "mapping process\n";
	
	my $mappingDir = "$dOut/mapping";
	mkdir $mappingDir unless (-d $mappingDir) ;
	
	## split genotype 
	{
		my ($lgFile) = glob("$dOut/grouping/*.lg");
		my $genotypeDir = "$mappingDir/genotype";
		mkdir $genotypeDir unless (-d $genotypeDir) ;

		print "perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir\n";
		print $log "perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir\n";
		`perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir`;
	}
	
	## write shell 
	my @genotypeFiles = glob("$mappingDir/genotype/*.genotype");
	open (SH,">$mappingDir/$fKey.cMap.sh") or die $!;
	foreach my $genotypeFile (@genotypeFiles) {
		my ($lg) = basename($genotypeFile) =~/(\w+)\.genotype$/ ;
		print SH "perl $Bin/sgsMapSmoothCycle.pl -i $genotypeFile -k $lg -d $mappingDir/$lg -mode 1 -n $cycle_t&&\n";
	}
	close (SH) ;
	
	## qsub 
	#my $hostname = `hostname`;chomp $hostname;
	#if ($hostname =~/cluster/) {
		print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
		print $log "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
		`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh`;
	#}else{
	#	print "ssh cluster -Y qsub-sge.pl --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
	#	print $log "ssh cluster -Y qsub-sge.pl --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
	#	`ssh cluster -Y qsub-sge.pl --reqsub --independent $mappingDir/$fKey.cMap.sh`;
	#}
	
	## check wheather all map files generated
	my @result = glob("$mappingDir/*/result");
	warn "mapping process unfinished, please check" if (@result != @genotypeFiles) ;

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
}

# step 3 : evaluate map 
if ($step == 3) {
	print "evaluate maps\n";
	print $log "evaluate maps\n";
	
	my $eMapDir = "$dOut/visualization";
	mkdir "$eMapDir" unless (-d $eMapDir) ;

	## prepare data
	`rm -rf $eMapDir/Data` if (-d "$eMapDir/Data") ;
	mkdir "$eMapDir/Data" unless (-d "$eMapDir/Data") ;
	`ln -s $dOut/mapping/*/result/*  $eMapDir/Data`;

	print "perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir\n";
	print $log "perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir\n";
	`perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir`;

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
}

# step4: reestimate mapdist using regression algorithm 
if ($step == 4) {
	print "reestimate map distance using regression mathods\n";
	print $log "reestimate map distance using regression mathods\n";
	
	my $edistDir = "$dOut/regression";
	mkdir "$edistDir" unless (-d $edistDir) ;
	mkdir "$edistDir/Data" unless (-d "$edistDir/Data") ;

	exit unless (-d "$dOut/visualization/Data") ;
	my @map = glob("$dOut/visualization/Data/*.map");
	`ln -s $dOut/visualization/Data/*.pwd $edistDir/Data `;
	`ln -s $dOut/visualization/Data/*.loc $edistDir/Data `;

	### write shell 
	open (SH,">$edistDir/$fKey.regression.sh") or die $!;
	foreach my $map (@map) {
		my ($key) = basename($map) =~/(\S+)\.map$/;
		print SH "perl $Bin/rearrangeFinalMap.pl -m $map -p $edistDir/Data/$key.pwd -k $key -d $edistDir/Data &&\n";
	}
	close (SH) ;

	### qsub
#	my $hostname = `hostname`;chomp $hostname;
#	if ($hostname =~/cluster/) {
		print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh\n";
		print $log "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh\n";
		`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh`;
#	}else{
#		print "ssh cluster -Y qsub-sge.pl --reqsub --independent $edistDir/$fKey.regression.sh\n";
#		print $log "ssh cluster -Y qsub-sge.pl --reqsub --independent $edistDir/$fKey.regression.sh\n";
#		`ssh cluster -Y qsub-sge.pl --reqsub --independent $edistDir/$fKey.regression.sh`;
#	}

	### estimate map
	print "perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir\n";
	print $log "perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir\n";
	`perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir `;


	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
print $log "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
close ($log) ;
# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub printLogHead {#
	my $fh=shift;
	my $user=`whoami`;chomp $user;
	my $program="$Bin/$Script";
	my $time=GetTime();
	my $currentPath=`pwd`;chomp $currentPath;
	my $job="perl $program"."\t".join("\t",@Original_ARGV)."\n";
	print $fh ";==========================================\n";
	print $fh join("\n",(";user: $user",";workDir : $currentPath",";job: $job")),"\n";
	print $fh ";==========================================\n\n";
}

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
  -i		<file>	Genotype file or loc file, forced
  -k		<str>	Key of output file,forced
  -d		<dir>	Directory where output file produced,optional,default [./]
  
  -n		<int>	The number of chromosomes,forced

  -minGroup	<int>	Minimum linkage group, optional, default [20]
  -maxGroup	<int>	Maximum linkage group, optional, default [500]
  
  -cycle	<int>	Maximum number of sgsMap&smooth cycle, optional, default [5]

  -step		<int>	Pipeline step
			<1>  determine linakge group 
			<2>  mapping process 
			<3>  evaluate maps 
			<4>  re-estimate map distance using regression method
  -only			Do only one pipeline step 
  -h		Help

USAGE
	print $usage;
	exit;
}