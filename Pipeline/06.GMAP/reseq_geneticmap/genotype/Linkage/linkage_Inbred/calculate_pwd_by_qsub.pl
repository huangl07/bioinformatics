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
my ($locfile,$dOut,$fKey,$onlyFirst,$minMarkerNum,$step,$Only,$maxproc,$interval,$resource,$type);
GetOptions(
				"help|?" =>\&USAGE,
				
				"i:s"=>\$locfile,
				"d:s"=>\$dOut,
				"k:s"=>\$fKey,
					
				"p:i"=>\$onlyFirst,
				"m:i"=>\$minMarkerNum,	
				
				"maxproc:i"=>\$maxproc,
				"interval:i"=>\$interval,
				"resource:s"=>\$resource,
				"type:s"=>\$type,

				"step:i"=>\$step,
				"only"=>\$Only,

				) or &USAGE;
&USAGE unless ($locfile and $fKey and $type );
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

my $file_type = `grep ^nloc $locfile` ? 'loc': 'genotype';
print "$file_type\n";
my $MARKER_NUM = $file_type eq 'genotype' ?  ( split(/\s+/,`less -S $locfile|grep ^$ -P -v |wc -l`))[0] : `grep ^nloc $locfile `;

$MARKER_NUM = ($MARKER_NUM =~ /(\d+)/)[0] ;
chomp $MARKER_NUM ;

#print $MARKER_NUM ;die;
if ($MARKER_NUM <= $minMarkerNum) {

	print "Total Marker Number is less than $minMarkerNum,calculate directly\n";

	print " perl $Bin/Calculation_Recombinant_Ratio_all.pl -i $locfile -o $dOut/$fKey.pwd -t $type\n";
	` perl $Bin/Calculation_Recombinant_Ratio_all.pl -i $locfile -o $dOut/$fKey.pwd -t $type `;

}else{

	print "Calcualte pairWise data via qsub\n";
	
	my $qsubDir = "$dOut/calculatePWDdir";
	mkdir $qsubDir unless (-d $qsubDir) ;

	my (%Data,%gFiles) = ();
	
	## read genotype file 



	########sort marker

	my %info;
	open (IN,"$locfile") or die $!;
	open (OUT,">$dOut/calculatePWDdir/sort_marker.loc") or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ or /^;/) ;
		if (/=/) {
			print OUT $_,"\n";
			next;
		}
		my ($marker,undef) = split ;
		$info{$marker} = $_ ;
	}

	map {print OUT $info{$_},"\n" } (sort {$a cmp $b} keys %info);

	close (IN) ;
	close (OUT) ;

################




	open (IN,"$dOut/calculatePWDdir/sort_marker.loc") or die $!;
	my $head ;
	my $order = 0 ;
	while (<IN>) {
		chomp;
		next if (/^$/ or /^;/) ;
		if (/=/) {
			$head .= "$_\n";
			next;
		}
		my ($marker) = split ;

		$Data{$marker}{'line'} = $_;
		$Data{$marker}{'order'} = ++$order;

	}
	close (IN) ;

	### split genotype file ###

	if ($step == 1) {

		print "split genotype file\n" ;

		my $partNum = int($MARKER_NUM/$onlyFirst + 1);

		for(1..$partNum){

			my $gFile = "$qsubDir/$fKey.part$_.loc";

			open (OUT,">$gFile") or die $!;
			print OUT $head;
			close (OUT) ;

			$gFiles{$_} = AbsolutePath('file',$gFile);

		}

		foreach my $marker (sort {$Data{$a}{'order'} <=> $Data{$b}{'order'}} keys %Data) {

			for(1..$partNum){

				if ($Data{$marker}{'order'} > ($_-1)*$onlyFirst) {

					my $gFile = "$qsubDir/$fKey.part$_.loc";

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
			mkdir ($dirname) if (!-d $dirname );
			if ($part == $END_PART) {

				print SH "perl $Bin/Calculation_Recombinant_Ratio_all.pl -i $gFiles{$part} -o $dirname/$basename.pwd -t $type &&\n";	

			}else{

				print SH "perl $Bin/Calculation_Recombinant_Ratio_all.pl -i $gFiles{$part} -o $dirname/$basename.pwd -t $type  -p $onlyFirst &&\n";
			
			}
		}

		close (SH) ;

		### qsub ###

		print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh \n";
		`sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh --maxproc $maxproc --interval $interval --resource $resource --reqsub $qsubDir/cpw.sh`;
		
		$step++ unless($Only);
		print "\n";
	}

	if ($step == 3) {

		print "merge all pwd files\n";

		` cat $qsubDir/*.pwd |grep '^;' -P -v  > $dOut/$fKey.pwd `;

	}

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

sub AbsolutePath
{		#..............
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

			This program calls sub progame : $Bin/Calculation_Recombinant_Ratio_all.pl
Usage:
  Options:
  -i		<file>	loc file, forced
  -d		<str>	Directory where output file produced,optional,default [./]
  -k		<str>	Key of output file,forced
  
  -p		<int>	integer skip number marker at split file, default 200
  -m		<int>	Minimum threshold of markers\' number for qsub, optional, default [400]
  -type             Ri1234...  

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



