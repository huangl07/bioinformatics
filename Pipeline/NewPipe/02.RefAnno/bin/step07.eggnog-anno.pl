#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$out,$type,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"type:s"=>\$type,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	) or &USAGE;
&USAGE unless ($fa and $out and $dsh );
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
 $type||="nuc";
my $emapper="emapper.py --translate --dmnd_db /mnt/ilustre/users/dna/Environment/biotools/eggnog-mapper/data/NOG.fa.dmnd --dbtype seqdb";
if ($type eq "proc") {
	$emapper="emapper.py";
}
open SH,">$dsh/step06.EGGNOG1.sh";
open In,$fa;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my $fname=basename($_);
	print SH $emapper." -m diamond  --data_dir /mnt/ilustre/users/dna/Environment/biotools/eggnog-mapper/data/ -o $out/$fname.eggnog -i $_\n";
}
close SH;
close In;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step06.EGGNOG1.sh --maxjob=20 --CPU 8";
`$job`;
open SH,">$dsh/step06.EGGNOG2.sh";
print SH "perl $Bin/bin/EGGanno.pl -i $fa -d $out/ -o $out/EGGNOG.anno\n";
close SH;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step06.EGGNOG2.sh";
`$job`;

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
  -ref	<file>	input genome name,fasta format,
  -gff	<file>	input genome gff file,
  -out	<dir>	output data prefix
  -chr	<file>	chromosome change file
  -dsh	<dir>	output work sh dir

  -h         Help

USAGE
        print $usage;
        exit;
}
