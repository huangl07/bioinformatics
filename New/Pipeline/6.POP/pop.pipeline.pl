#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro,$step,$chr,$bayes);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"gro:s"=>\$gro,
	"out:s"=>\$out,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
	"bayes:s"=>\$bayes,
	"step:s"=>\$step,
			) or &USAGE;
&USAGE unless ($vcf and $out);
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$step||=1;
$maf||=0.05;
$mis||=0.3;
$dep||=2;
open Log,">$out/work_sh/pop.$BEGIN_TIME.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "variant-filer \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.vcf-filter.pl -vcf $vcf -out $out/step01.vcf-filter -dsh $out/work_sh -maf $maf -mis $mis -dep $dep";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-filter Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==2) {
	print Log "########################################\n";
	print Log "structure-generic \n",my $time=time();
	print Log "########################################\n";
	$vcf=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.recode.vcf");
	my $job="perl $Bin/bin/step02.structure-generic.pl -pop $vcf -out $out/step02.structure-generic -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "structure-generic Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==3) {
	print Log "########################################\n";
	print Log "pca-generic \n",my $time=time();
	print Log "########################################\n";
	$vcf=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.recode.vcf");
	if ($gro) {
		$gro=ABSOLUTE_DIR("$gro");
	}else{
		$gro=ABSOLUTE_DIR("$out/step02.structure-generic/group.list");
	}
	my $job="perl $Bin/bin/step03.pca-generic.pl -vcf $vcf -out $out/step03.pca-generic -dsh $out/work_sh ";
	$job.="-gro $gro\n" if ($gro);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "pca-generic Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==4) {
	print Log "########################################\n";
	print Log "tree-generic \n",my $time=time();
	print Log "########################################\n";
	$vcf=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.recode.vcf");
	if ($gro) {
		$gro=ABSOLUTE_DIR("$gro");
	}else{
		$gro=ABSOLUTE_DIR("$out/step02.structure-generic/group.list");
	}
	my $job="perl $Bin/bin/step04.tree-generic.pl -vcf $vcf -out $out/step04.tree-generic -dsh $out/work_sh ";
	$job.="-gro $gro\n" if ($gro);
	$job.="-bayes " if ($bayes);
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "structure-generic Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==5) {
	print Log "########################################\n";
	print Log "kinship-generic \n",my $time=time();
	print Log "########################################\n";
	$vcf=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.recode.vcf");
	my $job="perl $Bin/bin/step05.kinship-generic.pl -vcf $vcf -out $out/step05.kinship-generic -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "structure-generic Done and elapsed time : ",time()-$time,"s\n";
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
	tree	pca		structure	kinship
	eg:
	perl $Script -i -o -k -c

Usage:

	-vcf	<file>	input vcf files
	-out	<dir>	output dir
	-gro	<file>	input group list
	-maf	<num>	maf default 0.05
	-mis	<num>	mis default 0.3
	-dep	<num>	default 2 
	-bayes		do bayes tree,default not
	-step		pipiline control
	-h			Help

USAGE
        print $usage;
        exit;
}
