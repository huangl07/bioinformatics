#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($est,$tpl,$out,$obs,$dsh,$number);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"est:s"=>\$est,
    "tpl:s"=>\$tpl,
    "obs:s"=>\$obs,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
    "number:s"=>\$number,
			) or &USAGE;
&USAGE unless ($est and $tpl and $obs and $dsh and $out);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$number||=30;

open SH,">$dsh/step02.fastsimcoal2.sh";
for(my $i=0;$i<=$number;$i++){
    my $workdir="$out/run$i";
    `mkdir $workdir`;
    `cp $obs/*.obs $workdir/`;
    `cp $tpl $workdir/`;
    `cp $est $workdir/`;
    print SH "cd $workdir/ && fastsimcoal2 -t pop.tpl -n100000 -N100000 -m -e pop.est -M 1e-5 -w 1e-5 -l 10 -L 40 -c 0 -q\n";
    #fastsimcoal2 -t $tpl -n100000 -N100000 -d -e $est -M 0.001 -l 10 -L 40 -q --multiSFS -C10
    }
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Maxjob 30 $dsh/step02.fastsimcoal2.sh";
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
	run fastsimcoal2
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -est	<file>	input est file 
  -tpl  <file>  input tpl file
  -obs  <dir>   input obs dir
  -number   <number>    defult 30
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
