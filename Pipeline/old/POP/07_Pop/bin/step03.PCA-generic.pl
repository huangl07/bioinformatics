#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($pop,$out,$dsh,$maf,$mis,$dep,$gro,$chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$pop,
	"gro:s"=>\$gro,
	"out:s"=>\$out,
	"chr:s"=>\$chr,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($pop and $out and $dsh and $chr );
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$pop=ABSOLUTE_DIR($pop);
open SH,">$dsh/step03.PCA-generic.sh";
open In,$pop;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	next if (!/pca/);
	my (undef,$id,$vcf)=split(/\s+/,$_);
	print SH "plink --vcf $vcf --pca 20 header tabs var-wts --out $out/$id.pca --allow-extra-chr --chr-set $chr && ";
	print SH "Rscript $Bin/bin/pca.R --infile $out/$id.pca.eigenvec --outfile $out/$id.pca --varfile $out/$id.pca.eigenval ";
	if ($gro) {
		print SH "--gro $gro\n"
	}else{
		print SH "\n";
	}
}
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step03.PCA-generic.sh";
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
  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
