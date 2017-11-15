#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$ref,$gff,$vcf,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"vcf:s"=>\$vcf,
	"dsh:s"=>\$dsh,
	"out:s"=>\$out,
	"proc:s"=>\$proc
			) or &USAGE;
&USAGE unless ($ref and $gff and $vcf and $out);
$dsh||=$out;
$proc||=20;
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
$vcf=ABSOLUTE_DIR($vcf);
$dsh=ABSOLUTE_DIR($dsh);
my $data="$out/data/ref";
open Config,">$out/snpEff.config";
print Config "data.dir =$out/data\nref.genome : ref\n";
close Config;
open SH,">$dsh/09-1.config.sh";
print SH "mkdir -p $out/data/ref/ &&";
print SH "ln -s $ref $out/data/ref/sequences.fa && ";
print SH "ln -s $gff $out/data/ref/genes.gff && ";
print SH "java -jar /mnt/ilustre/users/dna/.env//bin//snpEff.jar build -gff3 -v ref -c $out/snpEff.config";
close SH;
open SH,">$dsh/09-2.annovar.sh";
open ANNOLIST,">$out/anno.list";
open In,$vcf;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$vcfs)=split(/\s+/,$_);
	print SH "java -jar /mnt/ilustre/users/dna/.env//bin//snpEff.jar -v ref -csvStats $out/$id.anno.csv -c $out/snpEff.config $vcfs > $out/$id.anno.primary.vcf && vcftools --vcf $out/$id.anno.primary.vcf --recode-INFO ANN --recode --out $out/$id";
	print ANNOLIST "$id\t$out/$id.recode\t$out/$id.anno.csv\n";
}
close In;
close SH;
close ANNOLIST;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=30G -CPU 1 --maxjob $proc $dsh/09-1.config.sh";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --Resource mem=30G -CPU 1 --maxjob $proc $dsh/09-2.annovar.sh";
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
  -ref	<file>	input reference file 
  -gff	<file>	input gff file
  -proc <num>   number of process for qsub,default 20
  -vcf	<file>	input vcf list
  -out	<dir>	output dir
  -dsh	<dir>	output work_sh dir
  -h         Help

USAGE
        print $usage;
        exit;
}
