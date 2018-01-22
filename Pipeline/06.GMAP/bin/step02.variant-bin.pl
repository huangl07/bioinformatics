#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($marker,$dir,$dShell,$winsize,$stepsize,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"marker:s"=>\$marker,
				"popt:s"=>\$popt,
				"winsize:s"=>\$winsize,
				"stepsize:s"=>\$stepsize,
				"popt:s"=>\$popt,
				"out:s"=>\$dir,
				"dsh:s"=>\$dShell,
				) or &USAGE;
&USAGE unless ($marker and $dir and $dShell);
$marker=ABSOLUTE_DIR($marker);
mkdir $dir if (!-d $dir);
mkdir $dShell if (!-d $dShell);
$dir=ABSOLUTE_DIR($dir);
$dShell=ABSOLUTE_DIR($dShell);
$winsize||=500;
$stepsize||=$winsize/5;
$popt||="F2";
my $Key="Total";
my $step=1;
open Log,">$dShell/binmap.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "pesudo check \n",my $time=time();
	print Log "########################################\n";
	open SH,">$dShell/step02-1.pesudo.sh";
	if ($popt eq "CP") {
		print SH "perl $Bin/bin/binCP-pesudo.pl -i $marker -o $dir -k $Key";
	}else{
		print SH "perl $Bin/bin/binNOCP-pesudo.pl -i $marker -o $dir -k $Key";
	}
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dShell/step02-1.pesudo.sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "pesudo chr Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "wintype\n",my $time=time();
	print Log "########################################\n";
	open SH,">$dShell/step02-2.wintype.sh";
	print SH "perl $Bin/bin/binmap-wintype.pl -i $dir/$Key.male.matrix  -o $dir -k $Key.male -win $winsize -step $stepsize\n";
	print SH "perl $Bin/bin/binmap-wintype.pl -i $dir/$Key.female.matrix  -o $dir -k $Key.female -win $winsize -step $stepsize\n";
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dShell/step02-2.wintype.sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "wintype chr Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "phase\n",my $time=time();
	print Log "########################################\n";
	open SH,">$dShell/step02-3.phase.sh";
	print SH "perl $Bin/bin//binmap-phase.pl -i $dir/$Key.male.wintype -o $dir -k $Key.male\n";
	print SH "perl $Bin/bin//binmap-phase.pl -i $dir/$Key.female.wintype -o $dir -k $Key.female\n";
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $dShell/step02-3.phase.sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "phase Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "merge\n",my $time=time();
	print Log "########################################\n";
	open SH,">$dShell/step02-4.merge.sh";
	if ($popt eq "CP") {
		print SH "perl $Bin/bin/binCP-merge.pl -f $dir/$Key.female.bin.phase -m $dir/$Key.male.bin.phase -o $dir -k $Key";
	}else{
		print SH "perl $Bin/bin/binNOCP-merge.pl -f $dir/$Key.female.bin.phase -m $dir/$Key.male.bin.phase -o $dir -k $Key";
	}
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=3G --CPU 1 $dShell/step02-4.merge.sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "merge Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}


#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-marker	<file>	marker file
	-popt	<str>	popt type
	-out	<dir>	output file dir
	-dsh	<dir>	output shell dir
	-winsize	<num>	winsize [k]
	-stepsize	<num>	stepsize [k]
USAGE
	print $usage;
	exit;
}
