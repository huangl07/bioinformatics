#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$wpid,$mpid,$mbid,$wbid,$popt,$ann,$chr,$gff,);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"gff:s"=>\$gff,
	"popt:s"=>\$popt,
	"ref:s"=>\$ref,
	"out:s"=>\$out,
	"wp:s"=>\$wpid,
	"mp:s"=>\$mpid,
	"mb:s"=>\$mbid,
	"wb:s"=>\$wbid,
			) or &USAGE;
&USAGE unless ($vcf and $gff and $popt and $ref and $ann and $out and $mbid);
mkdir $out if (!-d $out);
$popt||="F2";
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$ann=ABSOLUTE_DIR($ann);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
$wpid||="";
$mpid||="";
$mbid||="";
my $dOut="$out/work_flow";
my $dSh="$out/work_sh";
my $step=1;
open LOG,">$dSh/Resequence.$BEGIN_TIME.log";
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "variant filter and format\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.variant-format.pl -vcf $vcf -out $dOut/01.filter -dsh $dSh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 2) {
	print LOG "########################################\n";
	print LOG "index and Gprime calc\n"; my $time=time();
	print LOG "########################################\n";
	my $variant="$dOut/pop.filtered.variant";
	my $job="perl $Bin/bin/step02.prime-calc.pl -variant $variant -out $dOut/02.prime -dsh $dSh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 3) {
	print LOG "########################################\n";
	print LOG "index calc\n"; my $time=time();
	print LOG "########################################\n";
	my $Iregion="$dOut/pop.index.region";
	my $Gregion="$dOut/pop.gprime.region";
	my $job="perl $Bin/bin/step03.gather-result.pl -iregion $Iregion -gregion $Gregion -out $dOut/03.result -dsh $outdir/work_sh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
close LOG;
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
  -vcf	<file>	input vcf file
  -out	<dir>	output dir
  -pid	<id>	input p id
  -bid	<id>	input b id
  -ann	<file>	input ann file
  -h         Help

USAGE
        print $usage;
        exit;
}
