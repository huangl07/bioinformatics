#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($smc,$group,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "smc:s"=>\$smc,
    "group:s"=>\$group,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($smc and $group and $dsh and $out);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step03.estimate.sh";

open IN,$group;
open OUT,">$out/smc.list";
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$pop)=split/\t/,$_;
    $hash{$pop}=1;
    }
close IN;
foreach my $pop (keys %hash){
    mkdir "$out/$pop" if (!-d "$out/$pop");
    print SH "/mnt/ilustre/users/dna/smcpp/build/bin/smc++ estimate -o $out/$pop  1.25e-8 $smc/$pop.chr*.smc.gz\n";
    print OUT "$pop\t$out/$pop/model.final.json\n";
    }

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Queue amd --Resource mem=10G --CPU 8  $dsh/step03.estimate.sh";
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
	run smc
	eg:
	perl $Script -group -chr -out -dsh -h

Usage:
  Options:
  -smc <dir>    the smc out file by step 2
  -group <file> the group list file 
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
