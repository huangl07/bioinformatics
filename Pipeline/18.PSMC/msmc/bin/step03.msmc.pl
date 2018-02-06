#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$out,$fg,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
    "g:s"=>\$fg,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($fin and $fg and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
    run msmc2
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<dir>	input msmc input file dir
  -g    <file>  group list file
  -out	<dir>  out dir 
  -dsh  <dir>   shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}

mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);

$fg=ABSOLUTE_DIR($fg);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step03.msmc.sh";
open OUT,">$out/msmc_output.list";
open IN,$fg;
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$pop)=split/\t/,$_;
    $hash{$pop}=1;
    }
close IN;


foreach my $pop (keys %hash){
    my @file=glob("$fin/*_$pop*");
    print SH "/mnt/ilustre/users/dna/.env/bin/msmc2 -t 16 -p 1*2+15*1+1*2  -I 0,1,2,3 -o $out/$pop\_output ",join(" ",@file),"\n";
    print OUT "$pop\t$out/$pop\_output\n";
    }

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=10G --CPU 16 $dsh/step03.msmc.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR
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

