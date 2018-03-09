#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($hapmap,$trait,$outdir);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"hapmap:s"=>\$hapmap,
	"trait:s"=>\$trait,
	"outdir:s"=>\$outdir
			) or &USAGE;
&USAGE unless ($hapmap and $trait and $outdir);
mkdir $outdir if (!-d $outdir);
mkdir "$outdir/work_sh" if (!-d "$outdir/work_sh");
open SH,">$outdir/work_sh/tassel-1.sh";
print SH "tassel -Xmx30G -fork1 -h $hapmap -KinshipPlugin -endPlugin -export $outdir/kinship\n";
print SH "tassel -Xmx30G -fork1 -h $hapmap -PrincipalComponentsPlugin -covariance true -endPlugin -export $outdir/pca\n";
close SH;
open SH,">$outdir/work_sh/tassel-2.sh";
print SH "tassel -fork1 -h $hapmap -fork2 -t $trait -fork3 -r $outdir/pca1.txt -excludeLastTrait -fork4 -k $outdir/kinship.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel Optimum -export $outdir/mlm\n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=30G --CPU 1 --maxjob 1 $outdir/work_sh/tassel-1.sh";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=30G --CPU 1 --maxjob 1 $outdir/work_sh/tassel-2.sh";
`$job`;

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
