#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$vcfanno,$anno,$out,$maf,$chr,$trt,$mis,$dep,$step,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
    "vcfanno:s"=>\$vcfanno,
    "anno:s"=>\$anno,
	"out:s"=>\$out,
    "trt:s"=>\$trt,
    "chr:s"=>\$chr,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
	"step:s"=>\$step,
    "stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $trt and $out);

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	gwas pipline
	eg:
	perl $Script -i -o -k -c

Usage:

-vcf    <file>	input vcf files
-vcfanno    <file> input vcf anno file
-anno   <file>  input pop anno file
-out    <dir>	output dir
-trt    <file>	the trit file
-maf    <num>	maf default 0.05
-mis    <num>	mis default 0.3
-dep    <num>	default 2
-chr    <num>   the chr number,count block by plink
-step   pipiline control
-stop   pipiline control
-h  Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$trt=ABSOLUTE_DIR($trt);
$vcfanno=ABSOLUTE_DIR($vcfanno);
$anno=ABSOLUTE_DIR($anno);
$step||=1;
$stop||=-1;
$maf||=0.05;
$mis||=0.3;
$dep||=2;

open Log,">$out/work_sh/pop.$BEGIN_TIME.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "vcf-filer and count block \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.vcf-filter.pl -vcf $vcf -chr $chr -trt $trt -out $out/step01.vcf-filter -dsh $out/work_sh -maf $maf -mis $mis -dep $dep";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "vcf-filter Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==2) {
	print Log "########################################\n";
	print Log "run GWAS \n",my $time=time();
	print Log "########################################\n";
	my $tritlist=ABSOLUTE_DIR("$out/step01.vcf-filter/trit/trit.list");
    my $tassellist=ABSOLUTE_DIR("$out/step01.vcf-filter/trit/trit.tassel.list");
    my $hmp=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.hapmap");
    my $job="perl $Bin/bin/step02.gwas.pl  -trt $tritlist -tassel $tassellist -hmp $hmp -out $out/step02.gwas/ -dsh $out/work_sh ";
    print Log "$job\n";
    `$job`;
    print Log "done!\n";
	print Log "########################################\n";
	print Log "GWAS Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
    $step++ if ($step ne $stop);
}
if ($step ==3) {
	print Log "########################################\n";
	print Log "gwas-result-combine \n",my $time=time();
	print Log "########################################\n";
    my $gwaslist=ABSOLUTE_DIR("$out/step02.gwas/gwas.dir.list");
    my $job="perl $Bin/bin/step03.gwas-combine.pl -list $gwaslist  -out $out/step03.gwas-combine -dsh $out/work_sh";
    print Log "$job\n";
    `$job`;
    print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "gwas-result-combine Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==4) {
	print Log "########################################\n";
	print Log "gwas-result \n",my $time=time();
	print Log "########################################\n";
    my $list=ABSOLUTE_DIR("$out/step03.gwas-combine/gwas.combine.list");
    my $block=ABSOLUTE_DIR("$out/step01.vcf-filter/plink.blocks");
    my $job="perl $Bin/bin/step04.gwas-result.pl -list $list -block $block -out $out/step04.gwas-result -dsh $out/work_sh ";
    print Log "$job\n";
   `$job`;
    print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "gwas-result Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==5) {
    print Log "########################################\n";
    print Log "gwas-region \n",my $time=time();
    print Log "########################################\n";
    my $list1=ABSOLUTE_DIR("$out/step04.gwas-result/Bonferroni1/Bonferroni1.list");
    my $list2=ABSOLUTE_DIR("$out/step04.gwas-result/Bonferroni2/Bonferroni2.list");
    my $job="perl $Bin/bin/step05.gwas-region.pl -list1 $list1 -list2 $list2 -anno $anno -vcf $vcfanno -out $out/step05.gwas-region -dsh $out/work_sh ";
    print Log "$job\n";
    `$job`;
    print Log "$job\tdone!\n";
    print Log "########################################\n";
    print Log "gwas-region Done and elapsed time : ",time()-$time,"s\n";
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

