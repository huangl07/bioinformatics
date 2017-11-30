#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use newPerlBase;
my $BEGIN_TIME=time();
my $Title="genetic_map:cp_linkage";
my $version="1.2.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fKey,$dOut,$minGroup,$maxGroup,$nchro,$cycle_t,$step,$only,$log,$ref,$modify,$sus,$diff_ratio);
my $test;
my $queue="middle.q";
my $cpu=50;
my $vf="5G";
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				
				"n:s"=>\$nchro,
				
				"ref:s"=>\$ref,
				"modify:s"=>\$modify,

				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,
				
				"diff_ratio:s"=>\$diff_ratio,
				"sus:s"=>\$sus,
				"cycle:s"=>\$cycle_t,

				"step:s"=>\$step,
				"only"=>\$only,
				) or &USAGE;
&USAGE unless ($fIn and $fKey and $nchro);
createLog($Title,$version,$$,"$dOut/log/",$test);	
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
$dOut = Cwd::abs_path($dOut);
$fIn = Cwd::abs_path($fIn);
$ref = Cwd::abs_path($ref) if (defined $ref and -f $ref);
$modify||=1;
$step||= 1;
$minGroup||=20;
$maxGroup||=500;
$cycle_t||=5;
$diff_ratio||=0.85;
$sus||=0.06;
open ($log,">","$dOut/$fKey.log") or die $!;
&printLogHead($log);

#-------------------------------------------------------------------
# Pipeline
#-------------------------------------------------------------------

# step 1 : determine linkage group 
if ($step == 1) {
	print "determine linkage groups\n";
	print $log "determine linkage groups\n";
	stepStart(1,"determine linkage groups");
	my $groupDir = "$dOut/grouping";

	my $job= "perl $Bin/grouping/linkage_group.pl -i $fIn -k $fKey -o $groupDir -n $nchro -minGroup $minGroup -maxGroup $maxGroup ";
	if (defined $ref && $modify == 1) {
		$job.=" -c $ref -reference -modify 1 ";
	}elsif (defined $ref) {
		$job.=" -c $ref";
	}
	runOrDie($job);
	print $log $job,"\n";
	print $job,"\n";
	#`$job`;
	

	die "determine linkage group failed" unless (glob("$dOut/grouping/*.lg")) ;

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
	stepTime(1);
}
if ($step == 2) {
	
	if (defined $ref && $modify == 1) {
		stepStart(2,"modify_by_reference");
		my ($lgFile) = glob("$dOut/grouping/*.lg");
		if ($modify == 1) {
			my $job="perl $Bin/modify/modify_by_reference.pl -i $ref -l $lgFile -g $fIn -p CP -o $dOut/modify -k $fKey";
			my $d_ratio = $diff_ratio-$diff_ratio*0.15;
			my $u_ratio	= $sus * 1.5;
			$job .= " -diff_ratio $d_ratio -sus $u_ratio ";
			print $log $job,"\n";
			print $job,"\n";
			runOrDie($job);
			#` $job `;
			$fIn=glob("$dOut/modify/*.final.genotype");
			stepTime(2);
		}
		$step++;
	}else{
		$step++;
	}
}



# step 3 : construct linkage map each linkage group 
if ($step == 3) {
	my $job;
	print "mapping process\n";
	print $log "mapping process\n";
	stepStart(3,"mapping process");
	my $mappingDir = "$dOut/mapping";
	mkdir $mappingDir unless (-d $mappingDir) ;
	
	## split genotype 
	{
		my ($lgFile) = glob("$dOut/grouping/*.lg");
		my $genotypeDir = "$mappingDir/genotype";
		mkdir $genotypeDir unless (-d $genotypeDir) ;
		print "perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir\n";
		print $log "perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir\n";
		$job="perl $Bin/splitGenotypeViaLG.pl -l $lgFile -g $fIn -d $genotypeDir";
		runOrDie($job);
	}

	## write shell 
	my @genotypeFiles = glob("$mappingDir/genotype/*.genotype"); #modify by qix 
	open (SH,">$mappingDir/$fKey.cMap.sh") or die $!;
	foreach my $genotypeFile (@genotypeFiles) {
		my ($lg) = basename($genotypeFile) =~/(\w+)\.genotype$/ ;
		print SH "perl $Bin/sgsMapSmoothCycle.pl -i $genotypeFile -k $lg -d $mappingDir/$lg -mode 1 -n $cycle_t -e $sus -diff_ratio $diff_ratio&&\n";
	}
	close (SH) ;
	## qsub 
	print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
	print $log "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh\n";
	#`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $mappingDir/$fKey.cMap.sh`;
	qsubOrDie("$mappingDir/$fKey.cMap.sh",$queue,$cpu,$vf);
	#&Check_qsub_error("$mappingDir/$fKey.cMap.sh");

	
	## check wheather all map files generated
	my @result = glob("$mappingDir/*/result");
	warn "mapping process unfinished, please check" if (@result != @genotypeFiles) ;

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
	stepTime(3);
}

# step 4 : evaluate map 
if ($step == 4) {
	my $job;
	print "evaluate maps\n";
	print $log "evaluate maps\n";
	stepStart(4,"evaluate maps");
	my $eMapDir = "$dOut/visualization";
	mkdir "$eMapDir" unless (-d $eMapDir) ;

	## prepare data
	`rm -rf $eMapDir/Data` if (-d "$eMapDir/Data") ;
	mkdir "$eMapDir/Data" unless (-d "$eMapDir/Data") ;
	`ln -s $dOut/mapping/*/result/*  $eMapDir/Data`;

	print "perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir\n";
	print $log "perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir\n";
	$job="perl $Bin/drawFlow/drawMapFlow.pl -id $eMapDir/Data -k $fKey -igeno $fIn -od $eMapDir";
	runOrDie($job);

	mkdir "$dOut/QTL_Data" if (!-d "$dOut/QTL_Data");
	mkdir "$dOut/QTL_Data/final/" if (!-d "$dOut/QTL_Data/final/");
	my @mapping_type=("male","female","sexAver");#modified by zhaogf
	for (my $i=0;$i<@mapping_type;$i++) {
		my @map=glob("$eMapDir/Data/*.$mapping_type[$i].map");
		open Out,">$dOut/QTL_Data/final/$fKey.$mapping_type[$i].final.map";
		my $markerNum=0;
		my %hash;
		foreach my $map (@map) {
			my $name=basename($map);
			$name=~s/\D+//g;
			print Out "group\t$name\n";
			open In,$map;
			while (<In>) {
				chomp;
				next if ($_ eq "" || /group/ || /^$/);
				print Out $_,"\n";			
				$markerNum++;
				my @line=split/\s+/,$_;
				$hash{$line[0]}=$line[1];
			}
			close In;
		}
		my @loc=glob("$eMapDir/Data/*.$mapping_type[$i].loc");
		open Out,">$dOut/QTL_Data/final/$fKey.$mapping_type[$i].final.loc";
		print Out "name = $fKey\n";
		print Out "popt = CP\n";
		print Out "nloc = $markerNum\n";
		my $t=0;
		foreach my $loc (@loc) {
			open In,$loc;
			while (<In>) {
				chomp;
				next if ($_ eq "" || /group/ || /^$/ || /^;/);
				if (/nind/) {
					print Out $_,"\n" if ($t == 0);
					$t=1;
				}else{
					next if (/\=/);
					my $marker=(split/\s+/,$_)[0];
					if (exists $hash{$marker}) {
						print Out $_,"\n";	
					}		
				}
			}
			close In;
		}
		close Out;
	}
   if ($ref && $modify == 1) {
	   my $job1="perl $Bin/drawAligmentRalationMap.pl -m $dOut/QTL_Data/rearrange/$fKey.male.final.map -a $ref -k $fKey -o $eMapDir/ ";
	   my $job2="perl $Bin/drawAligmentRalationMap.pl -m $dOut/QTL_Data/rearrange/$fKey.female.final.map -a $ref -k $fKey -o $eMapDir/ ";
	   my $job3="perl $Bin/drawAligmentRalationMap.pl -m $dOut/QTL_Data/rearrange/$fKey.sexAver.final.map -a $ref -k $fKey -o $eMapDir/ ";
	   runOrDie($job1);
	   runOrDie($job2);
	   runOrDie($job3);
	   print $log "$job1\n$job2\n$job3\n";
   }
	$step++ unless ($only) ;
	print "statistic linkage information\n";
	print "perl $Bin/stat_lg_info.pl -i $dOut/visualization/Data -o $dOut/visualization/ \n";
	print $log "perl $Bin/stat_lg_info.pl -i $dOut/visualization/Data -o $dOut/visualization/ \n";
	$job="perl $Bin/stat_lg_info.pl -i $dOut/visualization/Data -o $dOut/visualization/ "; #生成map.stat.xls modifed by qix
	runOrDie($job);
	print "\n";
	print $log "\n";
	stepTime(4);
}

# step5: reestimate mapdist using regression algorithm 
if ($step == 5) {
	stepStart(5,"reestimate map distance using regression mathods");
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
	print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh\n";
	print $log "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh\n";
#		`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --reqsub --independent $edistDir/$fKey.regression.sh`;
	### estimate map
	print "perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir\n";
	print $log "perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir\n";
#	`perl $Bin/drawFlow/drawMapFlow.pl -id $edistDir/Data -k $fKey -igeno $fIn -od $edistDir `;
	mkdir "$dOut/QTL_Data" if (!-d "$dOut/QTL_Data");
	mkdir "$dOut/QTL_Data/rearrange/" if (!-d "$dOut/QTL_Data/rearrange/");

	@map=glob("$edistDir/Data/*.map");
	open Out,">$dOut/QTL_Data/rearrange/$fKey.final.map";
	my $markerNum=0;
	foreach my $map (@map) {
		my $name=basename($map);
		next if($name=~/male/);   # added by shitw 2015/3/19 filter male and female map
		$name=(split(/\_/,$name))[0];
		$name=~s/\D+//g;
		print Out "group $name\n";
		open In,$map;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /group/ || /^$/);
			print Out $_,"\n";			
			$markerNum++;
		}
		close In;
	}
	my @loc=glob("$edistDir/Data/*.loc");
	open Out,">$dOut/QTL_Data/rearrange/$fKey.final.loc";
	print Out "name = $fKey\n";
	print Out "popt = CP\n";
	print Out "nloc = $markerNum\n";
	my $t=0;
	foreach my $loc (@loc) {
		my $name=basename($loc);
		next if($name=~/male/);   # added by shitw 2015/3/19 filter male and female map
		$name=~s/\D+//g;
		open In,$loc;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /group/ || /^$/ || /^;/);
			if (/nind/) {
				print Out $_,"\n" if ($t == 0);
				$t=1;
			}else{
				next if (/\=/);
				print Out $_,"\n";			
			}
		}
		close In;
	}
   if ($ref && $modify == 1) {
	   my $job="perl $Bin/drawAligmentRalationMap.pl -m $dOut/QTL_Data/rearrange/$fKey.final.map -a $ref -k $fKey.rearrange -o $edistDir/ ";
	   runOrDie($job);
	   print $log $job;
   }

	$step++ unless ($only) ;
	print "\n";
	print $log "\n";
	stepTime(5);
}
totalTime();	

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
sub Check_qsub_error {#
	# Check The qsub process if error happend 
	my $sh=shift;
	my @Check_file=glob "$sh*.qsub/*.Check";
	my @sh_file=glob "$sh*.qsub/*.sh";

	if ($#sh_file!=$#Check_file) {
		print "Their Some Error Happend in $sh qsub, Please Check..\n";
		die;
	}
	else {
		print "$sh qsub is Done!\n";
	}
}
sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription:
		diff_ratio 越低，纠错越粗暴
Usage:
  Options:
  -i		<file>	Genotype file or loc file, forced
  -k		<str>	Key of output file,forced
  -d		<dir>	Directory where output file produced,optional,default [./]
  -ref		<file>	input posi file
	-modify			modify genotype and linkage by reference default 1
  -n		<int>	The number of chromosomes,forced

  -diff_ratio	<float>	diff ratio,default 0.85;
  -sus		<float>	Average error ratio of genotype matrix,optional,default [0.06]

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
  -step     debug
USAGE
	print $usage;
	exit;
}
