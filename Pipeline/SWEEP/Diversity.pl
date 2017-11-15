#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$fgroup,$maf,$dep);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"group:s"=>\$fgroup,
	"maf:s"=>\$maf,
	"dep:s"=>\$dep,
	"out:s"=>\$dOut,
			) or &USAGE;
&USAGE unless ($vcf and $fgroup and $dOut);
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($vcf);
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "vcf filter\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step00.vcf-manage.pl -v $vcf -g $fgroup -out $dOut/00.vcf-manage -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "reference prepair\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 0) {
	print LOG "########################################\n";
	print LOG "parameter calc\n"; my $time=time();
	print LOG "########################################\n";
	$vcf="$dOut/00.vcf-manage/pop.filtered.recode.vcf";
	my $job="perl $Bin/bin/step01.parameter-calc.pl -vcf $vcf -out $dOut/01.paremeter-calc -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "reference prepair\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 0) {
	print LOG "########################################\n";
	print LOG "draw parameter\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step00.uniform.pl -fqlist $fqlist -out $dOut/00.uniform -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "reference prepair\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
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

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
