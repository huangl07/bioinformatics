#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$num,$order,$root,$maf,$mis,$dep,$group,$step,$chr,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"group:s"=>\$group,
    "num:s"=>\$num,
    "order:s"=>\$order,
    "root:s"=>\$root,
	"out:s"=>\$out,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
	"step:s"=>\$step,
    "stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($vcf and $order and $root and $out);
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
$num||="5";
$order=ABSOLUTE_DIR($order);
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
	print Log "run treemix \n",my $time=time();
	print Log "########################################\n";
	my $tre=ABSOLUTE_DIR("$out/step01.vcf-filter/pop.tmix.gz");
	my $job="perl $Bin/bin/step02.treemix.pl  -tre $tre  -num $num -root $root -out $out/step02.treemix -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "run treemix Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
    $step++ if ($step ne $stop);
}
if ($step ==3) {
	print Log "########################################\n";
	print Log "plot treemix \n",my $time=time();
	print Log "########################################\n";
	my $list=ABSOLUTE_DIR("$out/step02.treemix/treemix.list");
	my $job="perl $Bin/bin/step03.draw.pl -list $list -order $order -out $out/step03.plot -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "plot treemix Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	treemix pipline
	eg:
	perl $Script -i -o -k -c

Usage:

-vcf	<file>	input vcf files;must be give
-out	<dir>	output dir;must be give
-group	<file>	input group list;must be give
-num    <number>
-root   <string>
-order  <file>
-maf	<num>	maf default 0.05
-mis	<num>	mis default 0.3
-dep	<num>	default 2
-step   pipiline control
-stop   pipiline control
-h			Help

USAGE
        print $usage;
        exit;
}
