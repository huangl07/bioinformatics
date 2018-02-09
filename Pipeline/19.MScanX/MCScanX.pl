#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fgff1,$fIn2,$fgff2,$fOut,$fchr1,$fchr2,$step,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"p1|profa1:s"=>\$fIn1,
	"g1|refgff1:s"=>\$fgff1,
	"c1|chr1:s"=>\$fchr1,	# input chr.prefix
	"p2:s"=>\$fIn2,
	"g2:s"=>\$fgff2,
	"c2:s"=>\$fchr2,
	"o:s"=>\$fOut,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			); 		#;/or &USAGE;
&USAGE unless ($fIn1 and $fgff1 and $fchr1 and $fIn2 and $fgff2 and $fchr2 and $fOut );
$fIn1=&DIR($fIn1);
$fIn2=&DIR($fIn2);
$fgff1=&DIR($fgff1);
$fgff2=&DIR($fgff2);
$fchr1=&DIR($fchr1);
$fchr2=&DIR($fchr2);
mkdir $fOut if(!-d $fOut);
my $dir=&DIR($fOut);
mkdir "$dir/work_sh" if (!-d "$dir/work_sh");
$step||=1;
$stop||=-1;
# prepare
if($step==1){
	my $job1="perl $Bin/bin/1_1.blastp.pl -p1 $fIn1 -p2 $fIn2 -o $dir/1.Prepare -d $dir/work_sh ";
	`$job1`;
	
	my $job2="perl $Bin/bin/1_2.con_gff.pl -p1 $fIn1 -p2 $fIn2 -c1 $fchr1 -c2 $fchr2 -g1 $fgff1 -g2 $fgff2 -o -o $dir/1.Prepare -d -d $dir/work_sh";
	`$job2`;
	$step++ if ($step ne $stop);
}
# MCScanX p & produce .ctl
if($step==2){
	mkdir "$dir/2.MCScan" if (!-d "$dir/2.MCScan");
	my $blast=DIR("$dir/1.Prepare/p.blast");
	my $gff=DIR("$dir/1.Prepare/p.gff");
	`ln -s $dir/1.Prepare/p.blast $dir/2.MCScan/`;
	`ln -s $dir/1.Prepare/p.gff $dir/2.MCScan/`;	
	`perl $Bin/bin/2.png.pl -c1 $fchr1 -c2 $fchr2 -o $dir/2.MCScan -d $dir/work_sh`;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub DIR{
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;
		$dir=`pwd`;
		chomp $dir;
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
sub USAGE {
        my $usage=<<"USAGE";
Contact:        qingmei.cui\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -p1 -g1 -c1 -p2 -g2 -c2 -o
Usage:
  Options:
  -h         Help

USAGE
        print $usage;
        exit;
}
