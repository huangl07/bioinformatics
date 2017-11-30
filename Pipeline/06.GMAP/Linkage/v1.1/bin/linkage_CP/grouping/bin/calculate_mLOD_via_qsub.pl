#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="2.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($Genotype,$dOut,$fKey,$onlyFirst,$minMarkerNum,$step,$Only,$maxproc,$interval,$resource);
GetOptions(
				"help|?" =>\&USAGE,
				
				"i:s"=>\$Genotype,
				"d:s"=>\$dOut,
				"k:s"=>\$fKey,
					
				"p:i"=>\$onlyFirst,
				"m:i"=>\$minMarkerNum,	
				
				"maxproc:i"=>\$maxproc,
				"interval:i"=>\$interval,
				"resource:s"=>\$resource,

				"step:i"=>\$step,
				"only"=>\$Only,

				) or &USAGE;
&USAGE unless ($Genotype and $fKey);
#------------------------------------------------------------------
# Global parameters 
#------------------------------------------------------------------

$dOut||="./";

mkdir $dOut unless (-d $dOut) ;

$dOut = AbsolutePath('dir',$dOut);

$minMarkerNum||=400;
$onlyFirst||=200;

$maxproc||=50;
$interval||=2;
$resource||="vf=5G";

$step||=1;
# ------------------------------------------------------------------
# Pipeline
# ------------------------------------------------------------------

my $MARKER_NUM = (split(/\s+/,`less -S $Genotype|grep ^$ -P -v|wc -l `))[0]-1;

if ($MARKER_NUM <= $minMarkerNum) {

	print "Total Marker Number is less than $minMarkerNum,calculate directly\n";

	print "perl $Bin/calculateMLOD.pl -i $Genotype -o $dOut/$fKey.mLOD \n";

	` perl $Bin/calculateMLOD.pl -i $Genotype -o $dOut/$fKey.mLOD `;

}else{

	print "Calcualte pairWise data via qsub\n";
	
	my $qsubDir = "$dOut/calculate_mLOD_dir";
	mkdir $qsubDir unless (-d $qsubDir) ;

	my (%Data,%gFiles) = ();
	
	## read genotype file 
	
	open (IN,"$Genotype") or die $!;
	my $head = <IN>;
	my $order = 0;
	while (<IN>) {
		chomp;
		next if (/^$/ or /^;/) ;
		
		my ($marker) = split;

		$Data{$marker}{'line'} = $_;
		$Data{$marker}{'order'} = ++$order;
		
	}
	close (IN) ;

	### split genotype file ###

	if ($step == 1) {

		print "split genotype file\n";

		my $partNum = int($MARKER_NUM/$onlyFirst+1);

		for(1..$partNum){

			my $gFile = "$qsubDir/$fKey.part$_.genotype";

			open (OUT,">$gFile") or die $!;
			print OUT $head;	
			close (OUT) ;

			$gFiles{$_} = AbsolutePath('file',$gFile);

		}

		foreach my $marker (sort {$Data{$a}{'order'} <=> $Data{$b}{'order'}} keys %Data) {

			for(1..$partNum){

				if ($Data{$marker}{'order'} > ($_-1)*$onlyFirst) {

					my $gFile = "$qsubDir/$fKey.part$_.genotype";

					open (OUT,">>$gFile") or die $!;
					print OUT $Data{$marker}{'line'},"\n";		
					close (OUT) ;

				}
			}
		}

		$step++ unless ($Only) ;

		print "\n"

	}

	### write shell for qsub ###

	if ($step == 2) {

		print "Write shell for qsub\n";

		open (SH,">$qsubDir/cpw.sh") or die $!;
		my ($END_PART) = sort {$b <=> $a} keys %gFiles;

		foreach my $part (sort {$a <=> $b} keys %gFiles) {

			my $dirname = dirname($gFiles{$part});
			my $basename = basename($gFiles{$part});

			if ($part == $END_PART) {

				print SH "perl $Bin/calculateMLOD.pl -i $gFiles{$part} -o $dirname/$basename.mLOD &&\n";	

			}else{

				print SH "perl $Bin/calculateMLOD.pl -i $gFiles{$part} -o $dirname/$basename.mLOD -p $onlyFirst &&\n";
			
			}
		}

		close (SH) ;

		### qsub ###

#		my $hostname = `hostname`;
#		if ($hostname !~/cluster/) {

			print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh \n";
			`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh`;

#		}else{
#			print "qsub-sge.pl --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh \n";
#			`qsub-sge.pl --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh`;
#		}

		$step++ unless($Only);
		print "\n";
	}

	if ($step == 3) {

		print "merge all mLOD files\n";

		` cat $qsubDir/*.mLOD|grep '^Loci_1\tLoci_2\tMLOD' -P -v  > $dOut/body `;
		` cat $qsubDir/*.mLOD|grep '^Loci_1\tLoci_2\tMLOD' -P |head -1 > $dOut/headline `;
		` cat $dOut/headline $dOut/body > $dOut/$fKey.mLOD`;
		`rm -rf $dOut/body`;
		`rm -rf $dOut/headline`;
	
	}

}

`touch $dOut/calculate_mLOD.check`;
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

sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program:$0 
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Description: 
		
			Calculate pairWise Data via Qsub.

			This program calls sub progame : $Bin/calculateMLOD.pl
  Options:
  -i		<file>	Genotype file, forced
  -d		<str>	Directory where output file produced,optional,default [./]
  -k		<str>	Key of output file,forced
  
  -p		<int>	integer skip number marker at split file, default 200
  -m		<int>	Minimum threshold of markers\' number for qsub, optional, default [400]
  

  -maxproc	<int>	default 50
  -interval	<int>	default 10
  -resource	<str>	defualt "vf=5G"
  
  
  -only			pipeline step control.
  -step		<int>	pipeline beginning contol, default 1, means from beginning.
			
			1:Split Genotype file

			2:write to shell

			3:Merged ltcLod and uniq
  
  -h			Help

USAGE
	print $usage;
	exit;
}
