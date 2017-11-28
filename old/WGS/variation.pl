#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$dOut,$ref,$gff,$RAD,$step,$stop,$SV,$CNV,$realign);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"out:s"=>\$dOut,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"realign"=>\$realign,
	"RAD"=>\$RAD,
	"SV"=>\$SV,
	"CNV"=>\$CNV,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($fqlist and $dOut and $ref and $gff);
mkdir $dOut if (!-d $dOut);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
$fqlist=ABSOLUTE_DIR($fqlist);
$dOut=ABSOLUTE_DIR($dOut);
my ($bamlist,$chrlist,$gvcflist,$vcflist,$svlist,$cnvlist,$dict,$annolist,$metric,$dictfile);
mkdir "$dOut/work_sh" if (!-d "$dOut/work_sh");
$step||=0;
$stop||=-1;
open LOG,">$dOut/work_sh/Resequence.log";
my $prestep=$step;
if ($step == 0) {
	print LOG "########################################\n";
	print LOG "reference prepair\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step00.genome-config.pl -ref $ref -gff $gff -out $dOut/step00.genome-config -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	$dict=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.dict");
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	print LOG "########################################\n";
	print LOG "reference prepair\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "bwa mapping \n"; my $time=time();
	print LOG "########################################\n";
	$dict=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.dict");
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	my $job="perl $Bin/bin/step01.bwa-mapping.pl -ref $ref -fqlist $fqlist -out $dOut/step01.bwa-mapping -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "bwa mapping\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 2) {
	print LOG "########################################\n";
	print LOG "bam sort \n"; my $time=time();
	print LOG "########################################\n";
	$bamlist=ABSOLUTE_DIR("$dOut/step01.bwa-mapping/bam.list");
	my $job="perl $Bin/bin/step02.bam-sort.pl -bam $bamlist -out $dOut/step02.bam-sort -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "bam sort\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
	if ($RAD) {$step++}
}
if ($step == 3) {
	print LOG "########################################\n";
	print LOG "bwa mkdup \n"; my $time=time();
	print LOG "########################################\n";
	$bamlist=ABSOLUTE_DIR("$dOut/step02.bam-sort/bam.sort.list");
	my $job="perl $Bin/bin/step03.bam-mkdup.pl  -bam $bamlist -out $dOut/step03.bam-mkdup -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "bwa mkdup\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
	$step++ if (!$realign);
}
if ($step == 4) {
	print LOG "########################################\n";
	print LOG "bwa realign \n"; my $time=time();
	print LOG "########################################\n";
	if ($RAD) {
		$bamlist=ABSOLUTE_DIR("$dOut/step02.bam-sort/bam.sort.list");
	}else{
		$bamlist=ABSOLUTE_DIR("$dOut/step03.bam-mkdup/bam.mkdup.list");
	}
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	$dict=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.dict");
	my $job="perl $Bin/bin/step04.bam-realign.pl -bam $bamlist -out $dOut/step04.bam-realign -dsh $dOut/work_sh -ref $ref";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "bwa mkdup\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 5) {
	if (!$RAD) {
		$metric=ABSOLUTE_DIR("$dOut/step03.bam-mkdup/metric.list");
	}else{
		$metric="NULL";
	}
	$bamlist=ABSOLUTE_DIR("$dOut/step04.bam-realign/bam.realign.list");
	$chrlist=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.chrlist");
	$dictfile=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.dict");
	print LOG "########################################\n";
	print LOG "map stat \n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step12.mapstat.pl -bam $bamlist -chr $chrlist -out $dOut/step12.mapping-stat -dsh $dOut/work_sh -met $metric --dict $dictfile";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "map stat\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 6) {
	print LOG "########################################\n";
	print LOG "haplotyper \n"; my $time=time();
	print LOG "########################################\n";
	if ($realign) {
		$bamlist=ABSOLUTE_DIR("$dOut/step04.bam-realign/bam.realign.list");
	}else{
		if ($RAD) {
			$bamlist=ABSOLUTE_DIR("$dOut/step02.bam-sort/bam.sort.list");
		}else{
			$bamlist=ABSOLUTE_DIR("$dOut/step03.bam-mkdup/bam.mkdup.list");
		}
	}
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	$dict=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.dict");
	my $job="perl $Bin/bin/step05.haplotyper-nosplit.pl -bam $bamlist -ref $ref -dict $dict -out $dOut/step05.haplotyper -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "haplotyper\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 7) {
	print LOG "########################################\n";
	print LOG "gvcf combine\n"; my $time=time();
	print LOG "########################################\n";
	$gvcflist=ABSOLUTE_DIR("$dOut/step05.haplotyper/gvcf.list");
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	my $job="perl $Bin/bin/step06.gvcf-typing-nosplit.pl -gvcf $gvcflist -ref $ref -out $dOut/step06.gvcf-typing -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "gvcf combine\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
	$step++ if ($step ne $stop);
}
if ($step == 8) {
	print LOG "########################################\n";
	print LOG "gvcf typing \n"; my $time=time();
	print LOG "########################################\n";
	$gvcflist=ABSOLUTE_DIR("$dOut/step06.gvcf-typing/gvcf.combine.list");
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	my $job="perl $Bin/bin/step07.vcf-merge.pl -vcf $gvcflist -ref $ref -out $dOut/step07.vcf-merge -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "gvcf typing\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 9) {
	print LOG "########################################\n";
	print LOG "variant filter \n"; my $time=time();
	print LOG "########################################\n";
	$vcflist=ABSOLUTE_DIR("$dOut/step07.vcf-merge/pop.variant.vcf");
	$ref=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.fa");
	my $job="perl $Bin/bin/step08.vcf-filter.pl -vcf $vcflist -ref $ref -out $dOut/step08.vcf-filter -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "variant filter\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 10) {
	print LOG "########################################\n";
	print LOG "annovar \n"; my $time=time();
	print LOG "########################################\n";
	$vcflist=ABSOLUTE_DIR("$dOut/step08.vcf-filter/vcf.filter.list");
	my $dref=ABSOLUTE_DIR("$dOut/step00.genome-config/");
	my $job="perl $Bin/bin/step11.annovar-snpEff.pl -vcf $vcflist -dref $dref -out $dOut/step11.annovar -dsh $dOut/work_sh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "annovar\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ( $step == 11) {
	if ($SV) {
		print LOG "########################################\n";
		print LOG "sv calling \n"; my $time=time();
		print LOG "########################################\n";
		if ($realign) {
			$bamlist=ABSOLUTE_DIR("$dOut/step04.bam-realign/bam.realign.list");
		}else{
			if ($RAD) {
				$bamlist=ABSOLUTE_DIR("$dOut/step02.bam-sort/bam.sort.list");
			}else{
				$bamlist=ABSOLUTE_DIR("$dOut/step03.bam-mkdup/bam.mkdup.list");
			}
		}
		$gff=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.gff");
		my $job="perl $Bin/bin/step09.sv-calling.pl -bam $bamlist -out $dOut/step09.sv-calling -dsh $dOut/work_sh -gff $gff";
		print LOG "$job\n";
		`$job`;
		print LOG "$job\tdone!\n";
		$svlist=ABSOLUTE_DIR("$dOut/step09.sv-calling/sv.filter.list");
		print LOG "########################################\n";
		print LOG "sv calling\n Done and elapsed time : ",time()-$time,"s\n";
		print LOG "########################################\n";
	}
	$step++ if ($step ne $stop);
}
if ($step == 12) {
	if ($CNV) {
		print LOG "########################################\n";
		print LOG "cnv calling \n"; my $time=time();
		print LOG "########################################\n";
		if ($realign) {
			$bamlist=ABSOLUTE_DIR("$dOut/step04.bam-realign/bam.realign.list");
		}else{
			if ($RAD) {
				$bamlist=ABSOLUTE_DIR("$dOut/step02.bam-sort/bam.sort.list");
			}else{
				$bamlist=ABSOLUTE_DIR("$dOut/step03.bam-mkdup/bam.mkdup.list");
			}
		}
		$gff=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.gff");
		my $job="perl $Bin/bin/step10.cnv-calling.pl -bam $bamlist -out $dOut/step10.cnv-calling -dsh $dOut/work_sh -gff $gff";
		print LOG "$job\n";
		`$job`;
		print LOG "$job\tdone!\n";
		print LOG "########################################\n";
		print LOG "cnv calling\n Done and elapsed time : ",time()-$time,"s\n";
		print LOG "########################################\n";
	}
	$step++ if ($step ne $stop);
}
if ($step == 13) {
	$cnvlist=ABSOLUTE_DIR("$dOut/step10.cnv-calling/cnv.filter.list")if($CNV);
	$svlist=ABSOLUTE_DIR("$dOut/step09.sv-calling/sv.filter.list") if($SV);
	$vcflist=ABSOLUTE_DIR("$dOut/step08.vcf-filter/vcf.filter.list");
	$annolist=ABSOLUTE_DIR("$dOut/step11.annovar/anno.list");
	$gff=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.gff");
	$chrlist=ABSOLUTE_DIR("$dOut/step00.genome-config/ref.chrlist");
	print LOG "########################################\n";
	print LOG "variant stat \n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step13.variantstat.pl -vcf $vcflist -anno $annolist -out $dOut/step13.variant-stat -dsh $dOut/work_sh -gff $gff -chr $chrlist ";
	$job.=" -sv $svlist " if ($SV);
	$job.=" -cnv $cnvlist " if ($CNV);
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "variant stat \n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 14) {
	$cnvlist=ABSOLUTE_DIR("$dOut/step10.cnv-calling/cnv.filter.list")if($CNV);
	$svlist=ABSOLUTE_DIR("$dOut/step09.sv-calling/sv.filter.list") if($SV);
	$vcflist=ABSOLUTE_DIR("$dOut/step08.vcf-filter/vcf.filter.list");
	$annolist=ABSOLUTE_DIR("$dOut/step11.annovar/anno.list");
	print LOG "########################################\n";
	print LOG "report \n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step14.report.pl -vcf $vcflist  -mapstat $dOut/step12.mapping-stat -variant $dOut/step13.variant-stat -out $dOut/step14-report -anno $annolist";
	$job.=" -sv $svlist " if ($SV);
	$job.=" -cnv $cnvlist " if ($CNV);
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "report stat \n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
close LOG;
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
	eg:
	perl $Script -fqlist -out -ref -gff

Usage:
  Options:
  -fqlist	<file>	input file name
  -out	<dir>	output dir
  -ref	<file>	reference file
  -gff	<file>	gff file
  -sv	sv calling default off 
  -realign	indel realign default off
  -cnv	cnv calling default off 
  -RAD	RRL calling default off 
  -step	pipeline control
  -h         Help

USAGE
        print $usage;
        exit;
}
