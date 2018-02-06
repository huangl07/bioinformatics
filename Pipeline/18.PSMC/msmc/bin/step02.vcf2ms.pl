#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$out,$fg,$dsh,$missing);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
    "g:s"=>\$fg,
    "m:s"=>\$missing,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($fin and $fg and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
    build msmc format 
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input ms-hap name
  -g    <file>  group list file
  -m    <number> the misssing 
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
$fin=ABSOLUTE_DIR($fin);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open SH,">$dsh/step02.vcf2ms.sh";
open IN,$fg;
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($sample,$pop)=split/\t/,$_;
    push @{$hash{$pop}},$sample;
    }
close IN;


foreach my $pop (keys %hash){
    my $pops=join("\,",@{$hash{$pop}});
    print SH "cd $out && python $Bin/bin/make_input_MSMC_from_callsTab.py -i $fin -o $pop.msmc -m $missing -s \"$pops\" \n";
    }

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step02.vcf2ms.sh";
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

