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
my ($fqlist,$config,$dOut,$dShell);
GetOptions(
				"help|?" =>\&USAGE,
				"fqlist:s"=>\$fqlist,,
				"config:s"=>\$config,
				"out:s"=>\$dOut,
				"dsh:s"=>\$dShell,
				) or &USAGE;
&USAGE unless ($fqlist and $config and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$dOut=ABSOLUTE_DIR($dOut);
$dShell=ABSOLUTE_DIR($dShell);
$fqlist=ABSOLUTE_DIR($fqlist);
open In,$fqlist;
my %lane;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($laneID,$lib,$fq1,$fq2)=split(/\s+/,$_);
	push @{$lane{$laneID}{$lib}{fq1}},$fq1;
	push @{$lane{$laneID}{$lib}{fq2}},$fq2;
}
close In;
open In,$config;
open SH,">$dShell/step01.rawdata-prepair.sh";
open Out,">$dOut/fastq.list";
my %libID;
my %fq1;
my %fq2;
my %type;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	#LaneID ProjectID       LibID   LibType SampleID        SampleNeed      Enzyme  Barcode
	my ($runID,$laneID,$projectID,$libID,$libType,$majorID,$sampleID,$sampleNeed,$enzyme1,$enzyme2,$Barcode)=split(/\t/,$_);
	if (!exists $lane{$laneID}{$libID}) {
		die "error in fastqlist $laneID $libID not find";
	}
	for (my $i=0;$i<@{$lane{$laneID}{$libID}{fq1}};$i++) {
		my $fq1=$lane{$laneID}{$libID}{fq1}[$i];
		my $fq2=$lane{$laneID}{$libID}{fq2}[$i];
		if ($libType eq "PE") {
			print SH "ln -s $fq1 $dOut/$laneID:$libID:$projectID:$majorID:$sampleID.R1.fastq.gz\n";
			print SH "ln -s $fq2 $dOut/$laneID:$libID:$projectID:$majorID:$sampleID.R2.fastq.gz\n";
			print Out join("\t","$laneID:$libID:$projectID:$majorID:$sampleID","$dOut/$laneID:$libID:$projectID:$majorID:$sampleID.R1.fastq.gz","$dOut/$laneID:$libID:$projectID:$majorID:$sampleID.R2.fastq.gz",$enzyme1,$enzyme2),"\n";
		}else{
			my $info=join(":",$laneID,$libID,$projectID);
			$libID{$info}=$Barcode;
			$fq1{$info}=$fq1;
			$fq2{$info}=$fq2;
			$type{$info}=$libType;
			print Out join("\t","$laneID:$libID:$projectID:$majorID:$sampleID","$dOut/$laneID:$libID:$projectID\_$sampleID\_R1.fastq.gz","$dOut/$laneID:$libID:$projectID\_$sampleID\_R2.fastq.gz",$enzyme1,$enzyme2),"\n";
		}
	}
}
close Out;
foreach my $info (sort keys %libID) {
	if ($type{$info} eq "RAD") {
		print SH " /mnt/ilustre/users/dna/.env//bin/axe-demux  -z 6 -b $libID{$info} -t $dOut/$info.split.stat -f $fq1{$info} -r $fq2{$info} -F $dOut/$info -R $dOut/$info \n";
	}else{
		print SH " /mnt/ilustre/users/dna/.env//bin/axe-demux -z 6 -c -2 -b $libID{$info} -t $dOut/$info.split.stat -f $fq1{$info} -r $fq2{$info} -F $dOut/$info -R $dOut/$info \n";
	}
}
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dShell/step01.rawdata-prepair.sh  --Resource mem=7G --CPU 2  --Nodes 1 --Maxjob 19 ";
`$job`;
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
	-fqlist	<file>	input fastq list 
	-config	<file>	input config list
	-out	<dir>	output dir
	-dsh	<dir>	output shell dir
USAGE
	print $usage;
	exit;
}
