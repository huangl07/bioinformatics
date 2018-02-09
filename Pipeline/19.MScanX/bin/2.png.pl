#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fdsh,$dOut,$fchr1,$fchr2);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"c1:s"=>\$fchr1,
	"c2:s"=>\$fchr2,
	"o:s"=>\$dOut,	# first_class dir;
	"d:s"=>\$fdsh,
			); 		#;/or &USAGE;
&USAGE unless ($dOut and $fdsh);
mkdir $dOut if(!-d $dOut); # class one_dir
mkdir "$fdsh" if (!-d "$fdsh");
my $dir=DIR($dOut);
$fchr1=DIR($fchr1);
$fchr2=DIR($fchr2);
open SH3,">$fdsh/2.1.pctl.sh";
# MCScanX p && convert chrlists to 4 classes of ctls.
	print SH3 "cd $dir && ~dna/Environment/biotools/MCScanX/MCScanX/MCScanX p && perl $Bin/bin/ctl.pl -c1 $fchr1 -c2 $fchr2 -o $dir\n";
close SH3;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $fdsh/2.1.pctl.sh";
`$job`;
my $col=&DIR("$dir/p.html");
# java -cp /mnt/ilustre/users/dna/Environment/biotools/MCScanX/MCScanX/downstream_analyses circle_plotter 
open SH4,">$fdsh/2.2.png.sh";
	print SH4 "java -cp /mnt/ilustre/users/dna/Environment/biotools/MCScanX/MCScanX/downstream_analyses circle_plotter -g $dir/p.gff  -s $dir/p.collinearity -c $dir/pro.circle.ctl -o $dir/pro.circle.png\n";
	print SH4 "java -cp /mnt/ilustre/users/dna/Environment/biotools/MCScanX/MCScanX/downstream_analyses dot_plotter -g $dir/p.gff  -s $dir/p.collinearity -c $dir/pro.dot.ctl -o $dir/pro.dot.png\n";	
	print SH4 "java -cp /mnt/ilustre/users/dna/Environment/biotools/MCScanX/MCScanX/downstream_analyses dual_synteny_plotter -g $dir/p.gff  -s $dir/p.collinearity -c $dir/pro.dualsynteny.ctl -o $dir/pro.dualsynteny.png\n";	
	print SH4 "java -cp /mnt/ilustre/users/dna/Environment/biotools/MCScanX/MCScanX/downstream_analyses bar_plotter -g $dir/p.gff  -s $dir/p.collinearity -c $dir/pro.bar.ctl -o $dir/pro.bar.png\n";
close SH4;	
`perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $fdsh/2.2.png.sh`;
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
	perl $Script -c1 -c2 -o -d 
Usage:
  Options:
  -c1	<file>	input chrlist name
  -c2	<file>	output chrlist name
  -o	<dir>	output
  -d 	<work_sh>
  -h         Help

USAGE
        print $usage;
        exit;
}
