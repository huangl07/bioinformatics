#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"pop:s"=>\$vcf,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh );
$vcf=ABSOLUTE_DIR($vcf);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
my $sam="$out/sample.list";
open SH,">$dsh/step02.structure1.sh";
print SH "vcftools --vcf $vcf --plink --out $out/pop && ";
print SH "plink --file $out/pop --make-bed --out $out/pop --allow-extra-chr && ";
print SH "cut -f 1 -d \" \" $out/pop.fam  > $out/sample.list \n";
close SH;
open SH,">$dsh/step02.structure2.sh";
open List,">$out/structure.list";
for (my $i=2;$i<20;$i++) {
	print SH "cd $out && ";
	print SH "admixture $out/pop.bed $i --cv -j8 > $out/pop.$i.log && paste $sam $out/pop.$i.Q > $out/pop.$i.xls && ";
	print SH "Rscript $Bin/bin/structure.R --infile $out/pop.$i.xls --outfile $out/pop.$i\n";
	print List "$i\t$out/pop.$i.log\t$out/pop.$i.xls\n";
}
close SH;
close List;
open SH,">$dsh/step02.structure3.sh\n";
print SH "perl $Bin/bin/CVerror.pl -i $out/structure.list -o $out";
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --CPU 1 --Resource mem=3G  $dsh/step02.structure1.sh";
print $job,"\n";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --CPU 8 --Resource mem=20G $dsh/step02.structure2.sh";
print $job,"\n";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl --CPU 1 --Resource mem=3G  $dsh/step02.structure3.sh";
print $job,"\n";
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
  -pop	<file>	input pop list 
  -out	<dir>	output dir
  -dsh	<dir>	output work shell
  -h         Help

USAGE
        print $usage;
        exit;
}
