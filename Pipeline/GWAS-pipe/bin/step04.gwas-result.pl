#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list,$block,$out,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"list:s"=>\$list,
    "block:s"=>\$block,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($list and $block and $out and $dsh);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	retrive gwas result
	eg:
	perl $Script -i -b -o 

Usage:
  Options:
  -list	<file>	input  gwas combine list
  -block    <file>  input blocks file
  -out	<dir>
  -dsh  <dir>
  -h         Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
mkdir "$out/Bonferroni1" if (!-d "$out/Bonferroni1");
mkdir "$out/Bonferroni2" if (!-d "$out/Bonferroni2");
$list=ABSOLUTE_DIR($list);
$block=ABSOLUTE_DIR($block);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);


open IN,$list;
open OUT1,">$out/Bonferroni1/Bonferroni1.list";
open OUT2,">$out/Bonferroni2/Bonferroni2.list";
open SH,">$dsh/step04.gwas-result.sh";
my %region_list;

while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$gwas)=split/\t/,$_;
    mkdir "$out/Bonferroni1/" if (!-d "$out/Bonferroni1/");
    mkdir "$out/Bonferroni2/" if (!-d "$out/Bonferroni2/");
    print OUT1 "$trit\t$out/Bonferroni1/$trit.bf1.region\n";
    print OUT2 "$trit\t$out/Bonferroni2/$trit.bf2.region\n";
    print SH "perl $Bin/bin/gwas.result.pl -gwas $gwas -block $block -trit $trit -out1 $out/Bonferroni1/ -out2 $out/Bonferroni2/\n";
}

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.gwas-result.sh";
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

