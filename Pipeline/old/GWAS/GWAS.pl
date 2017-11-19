#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$trt,$maf,$mis,$dep,$gro,$step,$chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"trt:s"=>\$trt,
	"out:s"=>\$out,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"step:s"=>\$step,
			) or &USAGE;
&USAGE unless ($vcf and $out);
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$chr||=1;
$step||=0;
$maf||=0.05;
$mis||=0.3;
$dep||=2;
open Log,">$out/work_sh/GWAS.log";
if ($step == 0) {
	print Log "########################################\n";
	print Log "variant-filer \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step00.variant-filter.pl -vcf $vcf -out $out/step00.vcf-filter -dsh $out/work_sh -maf $maf -mis $mis -dep $dep";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-filter Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==1) {
	print Log "########################################\n";
	print Log "vcf-convert \n",my $time=time();
	print Log "########################################\n";
	$vcf=ABSOLUTE_DIR("$out/step00.vcf-filter/pop.filtered.recode.vcf");
	my $job="perl $Bin/bin/step01.variant-convert.pl -vcf $vcf -out $out/step01.vcf-convert -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "vcf-convert Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==2) {
	print Log "########################################\n";
	print Log "do GWAS \n",my $time=time();
	print Log "########################################\n";
	my $hap=ABSOLUTE_DIR("$out/step01.vcf-convert/pop.hapmap");
	my $job="perl $Bin/bin/step02.GWAS.pl -hap $hap -out $out/step02-GWAS -dsh $out/work_sh -trt $trt ";
	$job.="-gro $gro\n" if ($gro);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "do GWAS Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
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
  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -trt	<file>	input trt list
  -chr	<num>	scaffold num
  -maf	<num>	maf default 0.05
  -mis	<num>	mis default 0.3
  -dep	<num>	default 2 
  -step		pipiline control
  -h         Help

USAGE
        print $usage;
        exit;
}
