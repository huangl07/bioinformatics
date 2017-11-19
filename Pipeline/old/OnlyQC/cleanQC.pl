#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dOut,$fqlist,$projectID,$step);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"outdir:s"=>\$dOut,
	"fqlist:s"=>\$fqlist,
	"projectID:s"=>\$projectID,
	"step:s"=>\$step,
			) or &USAGE;
&USAGE unless ($dOut and $fqlist);
mkdir $dOut if (!-d $dOut);
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);
$step||=1;
mkdir "$dOut/work_sh" if (!-d "$dOut/work_sh");
mkdir "$dOut/01.CleanDataQC" if (!-d "$dOut/01.CleanDataQC");
mkdir "$dOut/02.QCreport" if (!-d "$dOut/02.QCreport");
my $tmp=time();
open Log,">$dOut/work_sh/QC.$tmp.log";

if ($step == 1) {
	print Log "########################################\n";
	print Log " check CleanData QC!"; my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.clean-qc.pl -fqlist $fqlist -out $dOut/01.CleanDataQC/ -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log " check CleanData QC Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "merge result!"; my $time=time();
	print Log "########################################\n";
	my $cleanlist=ABSOLUTE_DIR("$dOut/01.CleanDataQC/stat.list");
	my $cleanfig=ABSOLUTE_DIR("$dOut/01.CleanDataQC/fig/");
	my $job="perl $Bin/bin/step02.merge_report.pl  -cleanlist $cleanlist  -cleanfig $cleanfig  -out $dOut/02.QCreport/$$projectID.report.xls ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "merge result Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "summary result!"; my $time=time();
	print Log "########################################\n";
	my $fIn=ABSOLUTE_DIR("$dOut/02.QCreport/$projectID.qc-report.xls");
	my $job="perl $Bin/bin/step03.summary_report.pl -i $fIn  -o $dOut/02.QCreport/report.xls ";
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
  -fqlist	<dir>	input fastq list
  -outdir	<dir>	output file dir
  -projectID	<str>	projectID
  -step		<num>	step number
  -h         Help

USAGE
        print $usage;
        exit;
}
