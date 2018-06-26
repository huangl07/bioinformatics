#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list,$order,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"list:s"=>\$list,
    "order:s"=>\$order,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($list and $order);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$list=ABSOLUTE_DIR($list);
$order=ABSOLUTE_DIR($order);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step03.plot.sh";
open IN,$list;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($id,$file)=split/\t/,$_;
    print SH "cd $out && Rscript $Bin/bin/draw.treemix.r --infile $file --order $order --outfile $id\n";
    }
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step03.plot.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -fa  -out -dsh -h

Usage:
  Options:
  -list  <file>  input treemix out list
  -order <number>    the order file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
