#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "list:s"=>\$list,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($list and $dsh and $out);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step04.plot.sh";

open IN,$list;
my @smc;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($pop,$smc_out)=split/\t/,$_;
    push @smc,$smc_out;
    }
close IN;
my $pops=join(" ",@smc);
print SH "/mnt/ilustre/users/dna/smcpp/build/bin/smc++ plot $out/plot.png $pops\n";
print SH "/mnt/ilustre/users/dna/smcpp/build/bin/smc++ plot $out/plot.pdf $pops\n";

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.plot.sh";
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
	draw smc
	eg:
	perl $Script  -list -out -dsh -h

Usage:
  Options:
  -list <file>    the smc file list
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
