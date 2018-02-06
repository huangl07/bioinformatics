#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($psmc,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"psmcfa:s"=>\$psmc,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($psmc);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$psmc=ABSOLUTE_DIR($psmc);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step03.psmc.sh";
open OUT,">$out/psmc.list";
open IN,$psmc;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($pop,$psmcfa)=split/\t/,$_;
    print SH "/mnt/ilustre/users/dna/.env/bin/psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o $out/$pop.psmc $psmcfa && /mnt/ilustre/users/dna/.env/bin/psmc_plot.pl -p -R $out/$pop $out/$pop.psmc\n";
    print OUT "$pop\t$out/$pop.psmc\n";
    }
close OUT;
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step03.psmc.sh";
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
	perl $Script -psmcfa  -out -dsh -h

Usage:
  Options:
  -psmcfa  <file> list file  
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
