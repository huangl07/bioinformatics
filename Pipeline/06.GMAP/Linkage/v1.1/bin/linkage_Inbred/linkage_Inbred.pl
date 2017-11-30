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
my $Title="genetic_map:linkage_Inbred";
my $version="1.2.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($genotype,$dOut,$type,$ref,$modify,$fKey,$nChro,$record,$step,$MSTmap,$minGroup,$maxGroup,$bLOD,$eLOD,$sLOD,$firstN,$maxNloci ,$min_theh,$min_miss,$nCycle,$log);
my $test;
my $queue="middle.q";
my $cpu=50;
my $vf="5G";
GetOptions(
				"help|?" =>\&USAGE,
				## IO 
				"i:s"=>\$genotype,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut,
				## chro and population type 
				"nChro:s"=>\$nChro,
				"popt:s"=>\$type,
				
				"ref:s"=>\$ref,
				"modify:s"=>\$modify,

				## grouping 
				"minGroup:s"=>\$minGroup,
				"maxGroup:s"=>\$maxGroup,
				"b:s"=>\$bLOD,
				"e:s"=>\$eLOD,
				"s:s"=>\$sLOD,
				
				## mapping 
				"record"=>\$record,
				"MSTmap"=>\$MSTmap,
				
				"n:s"=>\$nCycle,
				## smooth
				"min_H:s"=>\$min_theh,
				"min_X:s"=>\$min_miss,

				"step:s"=>\$step,
				"test"=>\$test,
				) or &USAGE;
&USAGE unless ($genotype and $dOut and $nChro and $type and $fKey);
createLog($Title,$version,$$,"$dOut/log/",$test);	

mkdir ($dOut) if (!-d $dOut) ;
$genotype = Cwd::abs_path($genotype) ;
$dOut = Cwd::abs_path($dOut) ;
$ref = Cwd::abs_path($ref) if (defined $ref and -f $ref);


$modify||=1;

#==================================================================
#  Command line options 
#==================================================================

# var
$step ||= 1 ;
$nCycle||= 10 ;

# linkage Grouping
$bLOD ||= 3 ;
$eLOD ||= 20 ;
$sLOD ||= 1 ;
$minGroup ||= 20 ;
$maxGroup ||= 800 ;
$firstN ||= 100;
$maxNloci ||= 400;

# smooth para 
$min_theh  ||= 0.9 ;
$min_miss  ||= 0.8 ;

#==================================================================
# Mapping process 
#==================================================================

my $grouping_dir = "$dOut/grouping";
my $ordering_dir = "$dOut/mapping";
my $eMap_dir = "$dOut/visualization";

if ($step == 1) {
	print "determine linkage groups\n";
	stepStart(1,"determine linkage groups"); 
	my $groupDir = "$dOut/grouping";

	my $job= "perl $Bin/grouping/linkage_group.pl -i $genotype -k $fKey -o $groupDir -n $nChro -minGroup $minGroup -maxGroup $maxGroup ";
	if (defined $ref && $modify == 1) {
		$job.=" -c $ref -reference -modify 1 ";
	}elsif (defined $ref) {
		$job.=" -c $ref";
	}
	print $job,"\n";
	runOrDie($job);
	#`$job`;
	

	die "determine linkage group failed" unless (glob("$dOut/grouping/*.lg")) ;

	$step++;
	print "\n";
	stepTime(1);	
}
if ($step == 2) {
	if ($ref) {
		stepStart(2,"modify_by_reference"); 
		my ($lgFile) = glob("$dOut/grouping/*.lg");
		if ($modify == 1) {
			my $job="perl $Bin/modify/modify_by_reference.pl  -i $ref -l $lgFile -g $genotype -p $type -o $dOut/modify -k $fKey ";
			my $m_theh = $min_theh - $min_theh * 0.15;
			my $m_miss = $min_miss - $min_miss * 0.15;
			$job.=" -D $min_theh -M $min_miss ";
			print $job,"\n";
			runOrDie($job);
			#`$job `;
			$genotype=glob("$dOut/modify/*.final.genotype");
		}
		$step++;
		stepTime(2);
	}else{
		$step++;
	}
}


if ($step == 3) {

	print "-- marker ordering --\n";
	mkdir $ordering_dir unless (-d "$ordering_dir") ;
	stepStart(3,"marker ordering"); 
	my @lg = glob("$grouping_dir/*.lg");
	die "more than one linkage group file found in dir: $grouping_dir" if (@lg >= 2) ;
	
	my $job = "perl $Bin/split_genotype_by_LG.pl -i $genotype -l $lg[0] -d $ordering_dir/genotype -k $fKey -t $type";
	print "$job\n";
	runOrDie($job);
	#`$job`;

	#
	# wirte shell 
	#
	my @loc = glob("$ordering_dir/genotype/*.loc");
	open (SH,">$ordering_dir/$fKey.cMap.sh") or die $!;
	
	foreach my $gFile (@loc) {
		my ($lg) = $gFile =~/\.(\w+)\.loc$/;
		$job = "perl $Bin/map_smooth_cycle.pl -i $gFile -o $ordering_dir/$lg -k $lg -t $type -n $nCycle -min_H $min_theh -min_X $min_miss ";
		$job .= " -record" if ($record) ;
		$job .= " -MSTmap" if ($MSTmap || (!$record && !$MSTmap)) ; ## default MSTmap 

		########added by liangsh
		$job .= " -p $ref" if ($ref);

		print SH "$job &&\n";
	}

	close (SH) ;
	#
	# qsub 
	# 
#	my $hostname = `hostname`; chomp $hostname;
#	if ($hostname =~/cluster/) {
		qsubOrDie("$ordering_dir/$fKey.cMap.sh",$queue,$cpu,$vf);
		$job = "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh $ordering_dir/$fKey.cMap.sh --interval 20 --maxproc 50 --reqsub";
		print "$job\n";
		#`$job`;
#	}else{
#		$job = "ssh cluster -Y qsub-sge.pl $ordering_dir/$fKey.cMap.sh --interval 20 --maxproc 50 --reqsub";
#		print "$job\n";
#		`$job`;
#	}

	print "\n";
	$step++;
	stepTime(3);
}

if ($step == 4) {
	print "-- evaluate map --\n";
	mkdir $eMap_dir unless (-d "$eMap_dir") ;
	stepStart(4,"evaluate map"); 
	my %stastic ;
	$stastic{'record'} = "LG\tnloc\tnind\tsum\tmiss\tsingleton\tnoise\tmiss(%)\tsingleton(%)\tnoise(%)\tgap_5\tgap_5(%)\tdistance\taver.distance\n";
	$stastic{'MST'}    = "LG\tnloc\tnind\tsum\tmiss\tsingleton\tnoise\tmiss(%)\tsingleton(%)\tnoise(%)\tgap_5\tgap_5(%)\tdistance\taver.distance\n";

	my @lgDir =  glob("$ordering_dir/*/MSTmap");
	foreach my $lg (@lgDir) {
		my $lgBasename = basename($lg);
#		print "$lgBasename\n";$lgBasename=MST
		$lg=~s/\/MSTmap//g;#modified by zhaogf
		print "$lg\n";
		my $chr_ID=basename($lg);#modified by zhaogf
#		print "$chr_ID\n";
		mkdir "$eMap_dir/$lgBasename" unless (-d "$eMap_dir/$lgBasename") ;
		mkdir "$eMap_dir/$lgBasename/Data" unless (-d "$eMap_dir/$lgBasename/Data") ;
		mkdir "$eMap_dir/$lgBasename/Result" unless (-d "$eMap_dir/$lgBasename/Result") ;
		#
		# cp loc
		#
		#`cp $ordering_dir/genotype/*$lgBasename.loc $eMap_dir/$lgBasename/Data`;
		
		#
		# record 
		#
		if (-d "$lg/record") {
			my $xls = (glob "$lg/record/*xls")[0];
			my @temp;
			open (XLS, $xls) or die $!;
			while (<XLS>) {
				chomp;
				@temp = split ;
				next if ($temp[-1] !~ /\d+/);
#				$temp[0] = $lgBasename ;
				$temp[0] = $chr_ID ;#modified by zhaogf
				push @temp , $temp[-1] / $temp[1];
			}
			close (XLS) ;

			$stastic{'record'} .= join ("\t", @temp);
			$stastic{'record'} .= "\n";

			my @cycle = glob "$lg/record/cycle*" ;
			my $max_cycle = 0 ;

			foreach my $cycle ( @cycle ) {
				$max_cycle = $1 if ( $cycle=~ /cycle_(\d+)/ && $1 > $max_cycle) ;
			}

			` cp $lg/record/cycle_$max_cycle/*order.loc  $eMap_dir/$lgBasename/Data/ ` ;
			` cp $lg/record/cycle_$max_cycle/*.map  $eMap_dir/$lgBasename/Data/ ` ;
			` cp $lg/record/cycle_$max_cycle/*DH.pwd  $eMap_dir/$lgBasename/Data/` if ($type =~ /[Dd][Hh]/) ;
			` cp $lg/record/cycle_$max_cycle/*.pwd  $eMap_dir/$lgBasename/Data/` if ($type !~ /[Dd][Hh]/) ;
			` cp $lg/record/cycle_$max_cycle/map/*png $eMap_dir/$lgBasename/Result/ `;
		}
		#
		# MSTmap
		#
		if (-d "$lg/MSTmap/") {
			my $xls = (glob "$lg/MSTmap/*xls")[0];
			print "$xls\n";
			my @temp ;
			open ( XLS ,$xls ) or die $!;
			while (<XLS>) {
				chomp;
				@temp = split ;
				next if ($temp[-1] !~ /\d+/);
#				$temp[0] = $lgBasename;
				$temp[0] = $chr_ID ;
				push @temp , $temp[-1] / $temp[1];
			}
			close (XLS);

			$stastic{MST} .= join ("\t" , @temp );
			$stastic{MST} .= "\n";

			my @cycle = glob "$lg/MSTmap/cycle*";
			my $max_cycle = 0;
			foreach my $cycle ( @cycle ) {
				$max_cycle = $1 if ( $cycle=~ /cycle_(\d+)/ && $1 > $max_cycle) ;
			}
			`cp $lg/MSTmap/cycle_$max_cycle/*order.loc  $eMap_dir/$lgBasename/Data/  ` ;
			`cp $lg/MSTmap/cycle_$max_cycle/*.map  $eMap_dir/$lgBasename/Data/ ` ;
			`cp $lg/MSTmap/cycle_$max_cycle/*DH.pwd  $eMap_dir/$lgBasename/Data/`  if ($type =~ /[Dd][Hh]/ || $type =~ /BC\d+/i) ;  #马立祥
			`cp $lg/MSTmap/cycle_$max_cycle/map/*.png  $eMap_dir/$lgBasename/Result/ ` ;
			`cp $lg/MSTmap/cycle_$max_cycle/MST_ordering/*.mst.o $eMap_dir/$lgBasename/Data/`; #add by malx;
		}
	}

	if ($record){
		open (RES, ">$eMap_dir/$fKey.record.stat.xls") or die $! ;
		print RES $stastic{record} ;
		close (RES) ;
	}

	if ($MSTmap || (!$record && !$MSTmap)){
		open ( MTS , ">$eMap_dir/$fKey.MSTmap.stat.xls" ) or die $!  ;
		print MTS $stastic{MST} ;
		close (MTS) ;
	}
	
	mkdir "$dOut/QTL_Data/" if (!-d "$dOut/QTL_Data/");
	my @map=glob("$eMap_dir/*/Data/*MST*.map");
	open Out,">$dOut/QTL_Data/$fKey.final.map";
	my $markerNum=0;
	foreach my $map (@map) {
		my $name=basename($map);
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
	my @loc=glob("$eMap_dir/*/Data/*MST*.loc");
	open Out,">$dOut/QTL_Data/$fKey.final.loc";
	print Out "name = $fKey\n";
	print Out "popt = $type\n";
	print Out "nloc = $markerNum\n";
	my $t=0;
	foreach my $map (@loc) {
#		my $name=basename($map);
#		$name=~s/\D+//g;
		open In,$map;
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
	}close Out;
#----------------------------------------------------------
#modified by zhaogf
#	my @letter=("a".."z");
#	my @double_letter;
#	for(my $i=0;$i<@letter;$i++) {
#		for (my $j=0;$j<@letter ;$j++) {
#			my $double=join "",$letter[$i],$letter[$j];
#			push @double_letter,$double;
#		}
#	}
#modified by liangsh
open IN,"<$genotype" or die $!;
my @indi;
while (<IN>) {
	chomp;
	if (/Type|ID/) {
		(undef,undef,@indi) = split;
		last;
	}
}
close IN;
	open Out,">$dOut/QTL_Data/$fKey.R.final.loc";
	$t=0;
	foreach my $map (@loc) {
		open In,$map;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /group/ || /^$/ || /^;/);
			if (/nind/) {
				#my $indi=$_;
				#$indi=~s/\D+//g;
				my $head=join "\t",@indi;
				print Out "individual\t$head\n" if($t==0);
				$t=1;
			}else{
				next if (/\=/);
				print Out $_,"\n";			
			}
		}
		close In;
	}close Out;
#----------------------------------------------------------
   if ($ref && $modify == 1) {
	   my $job="perl $Bin/drawAligmentRalationMap.pl -m $dOut/QTL_Data/$fKey.final.map -a $ref -k $fKey -o $eMap_dir -png ";
	   runOrDie($job);
	   #`$job`;
	   print $job,"\n";
   }


	print "\n";
	$step++;
	stepTime(4);
}
totalTime();	
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

Usage:  近交群体连锁分析流程

  Options:

	-i		<file>		input genotype file, forced 
	-k		<str>		output file stem, forced 
	-d		<dir>		output directory, forced
	-ref		<file>	input posi file
		-modify			modify genotype and linkage by reference default 1

	-popt		<str>		population type, {BC1, DH, Ri\\d+}, forced
	-nChro		<str>		number of chromosome, forced
	-record				record mapping method
	-MSTmap				MSTmap mapping method, default 

	-minGroup	<int>		minimum group size allowed in linkage grouping, default [20]
	-maxGroup	<int>		maximum group size allowed in linkage grouping, default [800]
	
	-n		<int>		maximum cycles allowed, default [10]
	
	-min_H		<double>	smooth min threshold to change a loc, default [0.9]
	-min_X		<double>	smooth min misratio to change a loc, default [0.8]
	
	-step 				pipeline step, default [1]
					1 linkage grouping 
					2 marker ordering using MSTmap or record
					3 evaluate map
	-h				help
	-test			debug
USAGE

	print $usage;
	exit;

}
