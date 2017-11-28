#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro,$chr,$sam);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"chr:s"=>\$chr,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh and $chr);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
open In,$vcf;
open SH,">$dsh/step01.vcf-convert.sh";
open LIST,">$out/pop.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$vcf)=split(/\s+/,$_);
	print SH "vcftools --vcf $vcf --plink --out $out/$id && ";
	print SH "plink --file $out/$id --make-bed --out $out/$id --allow-extra-chr \n";
	print SH "perl $Bin/bin/sample.pl -i $vcf -o $out/sample.list\n";
	print SH "perl $Bin/bin/vcf2raxml.pl -i $vcf -o $out/$id.phylip\n";
	print LIST "tree\t$id\t$out/$id.phylip\n";
	print LIST "structure\t$id\t$out/$id.bed\n";
	print LIST "pca\t$id\t$vcf\n";
}
close In;
close LIST;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.vcf-convert.sh";
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
  -vcf	<file>	input vcf list
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -nchr	<num>	scaffold number 
  -h         Help

USAGE
        print $usage;
        exit;
}
