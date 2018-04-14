#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($trt,$vcf,$out,$chr,$dsh,$maf,$mis,$dep,$gro);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
    "chr:s"=>\$chr,
    "trt:s"=>\$trt,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
	"maf:s"=>\$maf,
	"mis:s"=>\$mis,
	"dep:s"=>\$dep,
			) or &USAGE;
&USAGE unless ($vcf and $chr and $trt and $out and $dsh );
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
mkdir "$out/trit" if (!-d "$out/trit");
$out=ABSOLUTE_DIR($out);
$trt=ABSOLUTE_DIR($trt);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);

$mis||=0.3;
$maf||=0.05;
$dep||=2;
$mis=1-$mis;
open SH,">$dsh/step01.vcf-filter.sh";
    print SH "/mnt/ilustre/users/dna/.env/bin/vcftools --remove-filtered-all --remove-indels --minDP $dep  --max-missing $mis --vcf $vcf --recode  --out $out/pop --maf $maf && $Bin/bin/vcf2hapmap.pl -i $out/pop.recode.vcf -o $out/pop.hapmap && cd $out && /mnt/ilustre/users/dna/.env/bin/plink --vcf $out/pop.recode.vcf  --blocks no-pheno-req --allow-extra-chr --chr-set $chr\n";
    print SH "$Bin/bin/trit.pl -trt $trt -out $out/trit/ \n";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step01.vcf-filter.sh";
`$job`;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
	filter snp and retrive fa 
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf files
  -chr  <number> the chr number
  -trt  <file>  input trit file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -maf	<num>	maf filter default 0.05
  -mis	<num>	mis filter default 0.3
  -dep	<num>	dep filter default 2

  -h         Help

USAGE
        print $usage;
        exit;
}
