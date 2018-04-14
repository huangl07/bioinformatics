#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list,$out,$dsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"list:s"=>\$list,
	"out:s"=>\$out,
    "dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($list  and $out and $dsh);
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
  -list	<file>	input  gwas list
  -out	<dir>
  -dsh  <dir>
  -h         Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);


$list=ABSOLUTE_DIR($list);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);


open IN,$list;

open SH,">$dsh/step03.gwas-combine.sh";
open OUT,">$out/gwas.combine.list";

while(<IN>){
    $_=~s/[\n\r]//g;
    my ($trit,$dir)=split/\t/,$_;
    print OUT "$trit\t$out/$trit.combine.txt\n";
    print SH "perl $Bin/bin/gwas.combine.pl -trit $trit -dir $dir -out $out\n";
}
close IN;
close OUT;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step03.gwas-combine.sh";
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

