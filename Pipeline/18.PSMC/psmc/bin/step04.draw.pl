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
	"list:s"=>\$psmc,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($psmc);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$psmc=ABSOLUTE_DIR($psmc);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step04.draw.sh";
open IN,$psmc;
my @id;
my @result;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($pop,$psmc)=split/\t/,$_;
    #print SH "/mnt/ilustre/users/dna/.env/bin/psmc_plot.pl -p -M pop1,pop3 pop-all pop1.psmc pop3.psmc\n";
    push(@id,$pop);
    push(@result,$psmc);
    }
print SH "/mnt/ilustre/users/dna/.env/bin/psmc_plot.pl -p -M ",join(",",@id)," $out/pop-all ",join(" ",@result),"\n";
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.draw.sh";
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
	perl $Script -list  -out -dsh -h

Usage:
  Options:
  -list  <file> list file  
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
