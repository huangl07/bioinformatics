#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($tre,$num,$root,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"tre:s"=>\$tre,
    "num:s"=>\$num,
    "root:s"=>\$root,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($tre and $root);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$tre=ABSOLUTE_DIR($tre);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$num||=5;
open SH,">$dsh/step02.treemix.sh";
open LI,">$out/treemix.list";
for(my $i=0;$i<=$num;$i++){
    print SH "/mnt/ilustre/users/dna/.env/bin/treemix -i $tre -m $i -root $root -o $out/pop_$i > $out/out.$i.log\n";
    print LI "pop_$i\t$out/pop_$i\n";
    }
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step02.treemix.sh";
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
  -tre  <file>  input treemix file
  -num  <number>    the number to run treemix;defoult run 5 times
  -root <string>  the treemix root parameter,comma-delimited list of populations to set on one side of the root (for migration)
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
