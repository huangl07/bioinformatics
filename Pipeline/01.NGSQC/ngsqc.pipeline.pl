#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fconfig,$dOut,$fqlist,$step,$stop,$runID);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"config:s"=>\$fconfig,
	"outdir:s"=>\$dOut,
	"fqlist:s"=>\$fqlist,
	"runID:s"=>\$runID,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($fconfig and $dOut and $fqlist);
mkdir $dOut if (!-d $dOut);
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);
mkdir "$dOut/work_sh" if (!-d "$dOut/work_sh");
$fconfig=ABSOLUTE_DIR($fconfig);
$step||=0;
$stop||=-1;
my $tmp=time();
open Log,">$dOut/work_sh/QC.$tmp.log";
if ($step == 0) {
	print Log "########################################\n";
	print Log "configure for all\n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step00.config.pl -config $fconfig -outdir $dOut/00.config -run $runID -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "configure Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ unless($step == $stop);
}
if ($step == 1) {
	print Log "########################################\n";
	print Log "split for RAD\n",my $time=time();
	print Log "########################################\n";
	my $config=ABSOLUTE_DIR("$dOut/00.config/$runID.config");
	my $job="perl $Bin/bin/step01.rawdata-prepair.pl -fqlist $fqlist -config $config -out $dOut/01.RawData/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "split Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ unless($step == $stop);
}
my %enzyme;
if ($step == 2) {
	print Log "########################################\n";
	print Log "NGSQC!"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/01.RawData/fastq.list");
	my $job="perl $Bin/bin/step02.rawdata-qc.pl -fqlist $fqlist -out $dOut/02.RawDataQC/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "NGSQC Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step==3){
	print Log "########################################\n";
	print Log " check RawData duplication!"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/01.RawData/fastq.list");
	my $job="perl $Bin/bin/step03.rawdata-dup.pl -fqlist $fqlist -out $dOut/03.RawDataDup/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "check RawData duplication Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step==4){
	print Log "########################################\n";
	print Log "fastq trim  !"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/01.RawData/fastq.list");
	my $job="perl $Bin/bin/step04.fastp-qc.pl -fqlist $fqlist -out $dOut/04.CleanData/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " fastq trim  Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step==5){
	print Log "########################################\n";
	print Log "NGSQC!"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/04.CleanData/fastq.list");
	my $job="perl $Bin/bin/step05.clean-qc.pl -fqlist $fqlist -out $dOut/05.CleanDataQC/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " NGSQC Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step==6){
	print Log "########################################\n";
	print Log " check CleanData duplication!"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/04.CleanData/fastq.list");
	my $job="perl $Bin/bin/step06.clean-dup.pl -fqlist $fqlist -out $dOut/06.CleanDataDup/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " check CleanData duplication Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step==7){
	#check pollution;
	print Log "########################################\n";
	print Log "check pollution!"; my $time=time();
	print Log "########################################\n";
	$fqlist=ABSOLUTE_DIR("$dOut/04.CleanData/fastq.list");
	my $job="perl $Bin/bin/Pollution/Pollution.pl -fqlist $fqlist -o $dOut/07.Pollution/ -k pollution -n 10000";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " Check pollution Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 8) {
	print Log "########################################\n";
	print Log "merge result!"; my $time=time();
	print Log "########################################\n";
	my $rawlist=ABSOLUTE_DIR("$dOut/02.RawDataQC/stat.list");
	my $rawdup=ABSOLUTE_DIR("$dOut/03.RawDataDup/dup.list");
	my $cleanlist=ABSOLUTE_DIR("$dOut/05.CleanDataQC/stat.list");
	my $cleanfig=ABSOLUTE_DIR("$dOut/05.CleanDataQC/fig/");
	my $cleandup=ABSOLUTE_DIR("$dOut/06.CleanDataDup/dup.list");
	my $pollution=ABSOLUTE_DIR("$dOut/07.Pollution/pollution.summary.xls");
	
	my $job="perl $Bin/bin/step08.merge_report.pl -rawlist $rawlist -rawdup $rawdup -cleandup $cleandup  -cleanlist $cleanlist -cleanfig $cleanfig  -pollution $pollution -o $dOut/08.QCreport/";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "merge result Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	
}
close Log;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
  -config	<file>	input file xls file name
  -fqlist	<dir>	input fastq list
  -outdir	<dir>	output file dir
  -runID	<str>	runID
  -step		<num>	step number
  -h         Help

USAGE
        print $usage;
        exit;
}
