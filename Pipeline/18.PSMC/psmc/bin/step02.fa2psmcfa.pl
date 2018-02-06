#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($fa);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$fa=ABSOLUTE_DIR($fa);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);

open SH,">$dsh/step02.fa2psmcfa.sh";
open OUT,">$out/psmcfa.list";
open IN,$fa;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($pop,$fasta)=split/\t/,$_;
    print SH "/mnt/ilustre/users/dna/.env/bin/fq2psmcfa -q20 $fasta  > $out/$pop.psmcfa\n";
    print OUT "$pop\t$out/$pop.psmcfa\n";
    }
close OUT;
close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step02.fa2psmcfa.sh";
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
	fa to psmcfa
	eg:
	perl $Script -fa  -out -dsh -h

Usage:
  Options:
  -fa  <file> fa list file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
