#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$out,$num,$chr,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"num:s"=>\$num,
	) or &USAGE;
&USAGE unless ($fa and $out and $dsh );
$num||=20;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir $dsh if (!-d $dsh);
$dsh=ABSOLUTE_DIR($dsh);
$fa=ABSOLUTE_DIR($fa);
open SH,">$dsh/step02.fasplit.sh";
print SH "perl $Bin/bin/splitFasta.pl -i $fa -o $out -n $num\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step02.fasplit.sh";
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
  -fa	<file>	input genome name,fasta format,
  -out	<dir>	output data prefix
  -dsh	<dir>	output work sh dir
  -num	<num>	input split number
  -h         Help

USAGE
        print $usage;
        exit;
}
