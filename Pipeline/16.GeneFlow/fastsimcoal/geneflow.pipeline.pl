#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$group,$step,$chr,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"group:s"=>\$group,
	"out:s"=>\$out,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
    "p:s"=>\$p,
	"step:s"=>\$step,
    "stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	population history by psmc
	eg:
	perl $Script -i -o -k -c

Usage:

	-vcf	<file>	input vcf files;must be give
	-out	<dir>	output dir;must be give
	-group	<file>	input group list;must be give
	-maf	<num>	maf default 0.05
	-mis	<num>	mis default 0.3
	-dep	<num>	default 2
    -p  default 4+25*2+4+6
    -step   pipiline control
    -stop   pipiline control
    -h			Help

USAGE
        print $usage;
        exit;
}
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$out=ABSOLUTE_DIR($out);
$vcf=ABSOLUTE_DIR($vcf);
$group=ABSOLUTE_DIR($group);
$step||=1;
$stop||=-1;
$maf||=0.05;
$mis||=0.3;
$dep||=2;

open Log,">$out/work_sh/pop.$BEGIN_TIME.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "variant-filer \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.vcf-filter.pl -vcf $vcf -group $group -out $out/step01.vcf-filter -dsh $out/work_sh -maf $maf -mis $mis -dep $dep";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "variant-filter Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==2) {
	print Log "########################################\n";
	print Log "fa-to-psmcfa \n",my $time=time();
	print Log "########################################\n";
	my $falist=ABSOLUTE_DIR("$out/step01.vcf-filter/fa.list");
	my $job="perl $Bin/bin/step02.fa2psmcfa.pl  -fa $falist  -out $out/step02.fa2psmcfa -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "vcf-to-psmcfa Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
    $step++ if ($step ne $stop);
}
if ($step ==3) {
	print Log "########################################\n";
	print Log "smc-estimate \n",my $time=time();
	print Log "########################################\n";
	my $psmcfa=ABSOLUTE_DIR("$out/step02.fa2psmcfa/psmcfa.list");
	my $job="perl $Bin/bin/step03.psmc.pl -psmcfa $psmcfa -p $p -out $out/step03.psmc -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "psmc Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==4) {
	print Log "########################################\n";
	print Log "draw-psmc \n",my $time=time();
	print Log "########################################\n";
	my $list=ABSOLUTE_DIR("$out/step03.psmc/psmc.list");
	my $job="perl $Bin/bin/step04.draw.pl -list $list -out $out/step04.draw -dsh $out/work_sh ";
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


