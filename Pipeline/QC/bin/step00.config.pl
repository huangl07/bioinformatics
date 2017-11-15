#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use utf8;
use Spreadsheet::ParseExcel;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fqdir,$dsh,$RunID);
GetOptions(
				"help|?" =>\&USAGE,
				"config:s"=>\$fIn,
				"outdir:s"=>\$dOut,
				"fqdir:s"=>\$fqdir,
				"runID:s"=>\$RunID,
				"dsh:s"=>\$dsh,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $dsh );
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step00.config-generic.sh";
print SH "perl $Bin/bin/readconfig.pl -config $fIn -outdir $dOut -run $RunID";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  $dsh/step00.config-generic.sh --Resource mem=3G";
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
	-config	<file>	input run xls file [forced]
	-outdir	<dir>	output dir,[forced]
	-fqdir	<dir>	"/mnt/ilustre/upload/hiseq/hiseq4000/20170816nXten/"
	-run	<num>	runID [forced]
	-dsh	<dir>	output work dir
USAGE
	print $usage;
	exit;
}
