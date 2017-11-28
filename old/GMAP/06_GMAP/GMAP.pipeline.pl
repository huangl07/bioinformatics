#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcflist,$annolist,$genelist,$step,$stop,$dOut,$PID,$BID,$ref,$hete,$chr);
my $runID;
my $fqlist;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcflist:s"=>\$vcflist,
	"nchr:s"=>\$nchr,
	"bin:s"=>\$bin,
	"PID:s"=>\$PID,
	"MID:s"=>\$MID,
	"Popt:s"=>\$Popt,
	"outdir:s"=>\$dOut,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcflist and $nchr and $PID and $MID and $outdir and $Popt);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir "$dOut/work_sh" if (!-d "$dOut/work_sh");
$vcflist=ABSOLUTE_DIR($vcflist);
my @PID=split(/\,/,$PID);
my @MID=split(/\,/,$MID);
$step||=0;
open Log,">$dOut/work_sh/BSA.log";
if ($step == 0) {
	print Log "########################################\n";
	print Log "variant-merge \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step00.variant-merge.pl -vcflist $vcflist -ref $ref -out $dOut/00.variant-merge -dsh $dOut/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 1) {
	print Log "########################################\n";
	print Log "bin map\n",my $time=time();
	print Log "########################################\n";
	my $job;
	if ($bin) {
		$job="perl $Bin/bin/step01.variant-bin.pl -marker $marker -position $pos -out $dOut/01.variant-bin -dsh $dOut/work_sh --winsize 500 --stepsize 100 -popt $Popt";
	}else{
		$job="mkdir $dOut/01.variant-bin/ && ln -s $marker $dOut/";
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "bin map Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}

if ($step == 3) {
	print Log "########################################\n";
	print Log "MLOD calculated\n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step02.linkage-group.pl -marker $marker -position $pos -out $dOut/02.linkage-group -dsh $dOut/work_sh\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "MLOD calculated Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "linkage group\n",my $time=time();
	print Log "########################################\n";
	#my $job="perl $Bin/bin/step01.index_calc.pl -vcf $vcf -PID $PID -BID $BID -out $dOut/01.index-calc/ -dsh $dOut/work_sh -anno $anno -gene $genelist";
	$job.=" -hete " if ($hete);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "linkage group Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 5) {
	print Log "########################################\n";
	print Log "linkage phase\n",my $time=time();
	print Log "########################################\n";
	#my $job="perl $Bin/bin/step01.index_calc.pl -vcf $vcf -PID $PID -BID $BID -out $dOut/01.index-calc/ -dsh $dOut/work_sh -anno $anno -gene $genelist";
	$job.=" -hete " if ($hete);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "linkage phase Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 6) {
	print Log "########################################\n";
	print Log "marker order\n",my $time=time();
	print Log "########################################\n";
	#my $job="perl $Bin/bin/step01.index_calc.pl -vcf $vcf -PID $PID -BID $BID -out $dOut/01.index-calc/ -dsh $dOut/work_sh -anno $anno -gene $genelist";
	$job.=" -hete " if ($hete);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "marker order Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 7) {
	print Log "########################################\n";
	print Log "map evolution\n",my $time=time();
	print Log "########################################\n";
	#my $job="perl $Bin/bin/step01.index_calc.pl -vcf $vcf -PID $PID -BID $BID -out $dOut/01.index-calc/ -dsh $dOut/work_sh -anno $anno -gene $genelist";
	$job.=" -hete " if ($hete);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "map evolution Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}


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
  -vcflist	<file>	input vcf list
  -annolist	<file>	input anno list
  -genelist	<file>	input gene summary file
  -ref	<file>	reference fasta file
  -chr	<file>	reference chr list
  -PID	<str>	parent ID
  -BID	<str>	bulk ID
  -hete		for F1 default off
  -outdir	<dir>	output file dir
  -step		<num>	step number
  -h         Help

USAGE
        print $usage;
        exit;
}
