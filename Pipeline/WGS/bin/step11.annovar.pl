#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$vcflist,$dOut,$dShell,$sv,$cnv,$dref,$gff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcflist,
	"dref:s"=>\$dref,
	"proc:s"=>\$proc,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell
			) or &USAGE;
&USAGE unless ($vcflist and $dref and $dOut and $dShell);
$proc||=20;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
$vcflist=ABSOLUTE_DIR($vcflist);
mkdir $dShell if (!-d $dShell);
$dShell=ABSOLUTE_DIR($dShell);
$dref=ABSOLUTE_DIR($dref);
open SH,">$dShell/step11.annovar1.sh";
open Out,">$dOut/anno.list";
open In,$vcflist;
my $samplelist;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($type,$vcf)=split(/\s+/,$_);
	die $vcf if(!-f $vcf);
	open Vcf,$vcf;
	while (<Vcf>) {
		chomp;
		next if ($_ eq "" || /^##/);
		if (/^#/) {
			my ($CHR,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@sample)=split(/\s+/);
			$samplelist=join(",",@sample);
		}
	}
	close Vcf;
	print SH "perl $Bin/bin/annovar/convert2annovar.pl -format vcf4 -filter pass --includeinfo $vcf -allsample --outfile $dOut/$type 2>$dOut/$type.convert.log \n";
	print SH "perl $Bin/bin/annovar/convert2annovar.pl -format vcf4 -filter pass --includeinfo $vcf -allsample --withfreq --outfile $dOut/$type.pop.avinput 2>$dOut/$type.convert.log \n";
}
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=3G --CPU 1 --maxjob $proc  $dShell/step11.annovar1.sh";
print $job;
`$job`;
open SH,">$dShell/step11.annovar2.sh";
my @avinput=glob("$dOut/*.avinput");
foreach my $avinput (@avinput) {
	my ($type,$sample,undef)=split(/\./,basename($avinput));
	print SH "perl $Bin/bin/annovar/annotate_variation.pl --buildver ref $avinput $dref --outfile $dOut/$type.$sample >$dOut/$type.$sample.anno.log 2>&1 &&";
	if ($sample eq "pop") {
		print SH "perl $Bin/bin/annovar/annotation_pop.pl -exonic $dOut/$type.$sample.exonic_variant_function -variant $dOut/$type.$sample.variant_function -outfile $dOut/$sample.$type.anno -indi $samplelist\n";
	}else{
		print SH "perl $Bin/bin/annovar/annotation.pl -exonic $dOut/$type.$sample.exonic_variant_function -variant $dOut/$type.$sample.variant_function -outfile $dOut/$sample.$type.anno -indi $sample\n";
	}
	print Out "$sample\t$type\t$dOut/$sample.$type.anno\n";
}
close SH;
close Out;
$job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=3G --CPU 1 --maxjob $proc $dShell/step11.annovar2.sh";
print $job;
`$job`;
open SH,">$dShell/step11.annovar3.sh";
print SH "perl $Bin/bin/mergeAnno.pl -i $dOut/anno.list -o $dOut/pop.merge.anno -l $samplelist \n";
close SH;
$job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl --Resource mem=3G --CPU 1 --maxjob $proc  $dShell/step11.annovar3.sh";
print $job;
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
  -vcf	<file>	input vcflist file
  -dref	<file>	input dref file
  -proc <num>	number of process for qsub,default 20
  -out	<dir>	output dir 
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
