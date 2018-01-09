#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$ref,$gff,$RAD,$step,$stop,$SV,$CNV,$realign,$anno);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"outdir:s"=>\$outdir,
	"anno:s"=>\$anno,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"RAD"=>\$RAD,
	"sv"=>\$SV,
	"cnv"=>\$CNV,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($fqlist and $outdir and $ref and $gff);
$anno||="./";
mkdir $outdir if (!-d $outdir);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
$anno=ABSOLUTE_DIR($anno);
$fqlist=ABSOLUTE_DIR($fqlist);
$outdir=ABSOLUTE_DIR($outdir);
mkdir "$outdir/work_sh" if (!-d "$outdir/work_sh");
$step||=1;
$stop||=-1;
open LOG,">$outdir/work_sh/Resequence.$BEGIN_TIME.log";
my $prestep=$step;
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "fastq qc\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastq-qc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 2) {
	print LOG "########################################\n";
	print LOG "reference prepair\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step02.ref-config.pl -ref $ref -gff $gff  -out $outdir/02.ref-config -dsh $outdir/work_sh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 3) {
	print LOG "########################################\n";
	print LOG "bwa mapping\n"; my $time=time();
	print LOG "########################################\n";
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step03.bwa-mapping.pl -ref $ref -fqlist $fqlist  -out $outdir/03.bwa-mapping -dsh $outdir/work_sh  -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 4) {
	print LOG "########################################\n";
	print LOG "bam sort\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/03.bwa-mapping/bam.list");
	my $job="perl $Bin/bin/step04.bam-sort.pl -bam $bamlist -out $outdir/04.bam-sort -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
	if ($RAD) {
		$step++;
	}
}
if ($step == 5) {
	print LOG "########################################\n";
	print LOG "bam mkdup\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/04.bam-sort/bam.list");
	my $job="perl $Bin/bin/step05.bam-mkdup.pl -bam $bamlist -out $outdir/05.bam-mkdup -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 6) {
	print LOG "########################################\n";
	print LOG "map-stat\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/05.bam-mkdup/bam.list");
	if ($RAD) {
		$bamlist=ABSOLUTE_DIR("$outdir/04.bam-sort/bam.list");
	}
	my $chr=ABSOLUTE_DIR("$outdir/02.ref-config/ref.chrlist");
	my $metric=ABSOLUTE_DIR("$outdir/05.bam-mkdup/metric.list");
	my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step06.map-stat.pl -bam $bamlist -met $metric -dict $dict -chr $chr -out $outdir/06.map-stat -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 7) {
	print LOG "########################################\n";
	print LOG "haplotype\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/05.bam-mkdup/bam.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step07.haplotyper.pl -bam $bamlist -ref $ref -dict $dict -out $outdir/07.haplotype -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}

if ($step == 8) {
	print LOG "########################################\n";
	print LOG "gvcf typing\n"; my $time=time();
	print LOG "########################################\n";
	my $gvcf=ABSOLUTE_DIR("$outdir/07.haplotype/gvcf.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step08.gvcf-typing.pl -dict $dict -gvcf $gvcf -ref $ref  -out $outdir/08.gvcf-typing -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 9) {
	print LOG "########################################\n";
	print LOG "vcf-filter\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/08.gvcf-typing/pop.variant.vcf");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step09.vcf-filter.pl -vcf $vcf -ref $ref -out $outdir/09.vcf-filter -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 10) {
	print LOG "########################################\n";
	print LOG "annovar\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/09.vcf-filter/vcf.list");
	my $con=ABSOLUTE_DIR("$outdir/02.ref-config/snpEff.config");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step10.annovar.pl -con $con -vcf $vcf -ref $ref -anno $anno -out $outdir/10.annovar -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
	if (!$SV) {
		$step++;
	}
	if (!$CNV) {
		$step++;
	}
}
if ($step == 11) {
	print LOG "########################################\n";
	print LOG "sv call\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/05.bam-mkdup/bam.list");
	my $gff=ABSOLUTE_DIR("$outdir/02.ref-config/ref/genes.gff");
	my $job="perl $Bin/bin/step11.sv-call.pl -bam $bamlist -gff $gff -out $outdir/11.sv-call -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 12) {
	print LOG "########################################\n";
	print LOG "cnv call\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/05.bam-mkdup/bam.list");
	my $job="perl $Bin/bin/step12.cnv-call.pl -bam $bamlist  -gff $gff -out $outdir/12.cnv-call -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 13) {
	print LOG "########################################\n";
	print LOG "variant stat\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/10.annovar/anno.list");
	my $chr=ABSOLUTE_DIR("$outdir/02.ref-config/ref.chrlist");
	my $gff=ABSOLUTE_DIR("$outdir/02.ref-config/ref/genes.gff");
	my $job="perl $Bin/bin/step13.var-stat.pl -vcf $vcf -chr $chr -gff $gff  -out $outdir/13.variant-stat -dsh $outdir/work_sh -proc 20 ";
	if ($SV) {
		my $sv=ABSOLUTE_DIR("$outdir/11.sv-call/sv.filter.list");
		$job.=" -sv $sv ";
	}
	if ($CNV) {
		my $cnv=ABSOLUTE_DIR("$outdir/12.cnv-call/cnv.filter.list");
		$job.=" -cnv $cnv ";
	}
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 14) {
	print LOG "########################################\n";
	print LOG "report generic\n"; my $time=time();
	print LOG "########################################\n";
	my $fqdir=ABSOLUTE_DIR("$outdir/01.fastq-qc");
	my $mapstat=ABSOLUTE_DIR("$outdir/06.map-stat");
	my $varstat=ABSOLUTE_DIR("$outdir/13.variant-stat");
	my $vcf=ABSOLUTE_DIR("$outdir/10.annovar/anno.list");
	my $anno=ABSOLUTE_DIR("$outdir/10.annovar");
	my $job="perl $Bin/bin/step14.report.pl -fastqc $fqdir -mapstat $mapstat -variant $varstat -vcf $vcf -out $outdir/14.report -annotate $anno";
	if ($SV) {
		my $sv=ABSOLUTE_DIR("$outdir/11.sv-call/sv.filter.list");
		$job.=" -sv $sv ";
	}
	if ($CNV) {
		my $cnv=ABSOLUTE_DIR("$outdir/12.cnv-call/cnv.filter.list");
		$job.=" -cnv $cnv ";
	}
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
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

	-fqlist	<file>	input file name
	-outdir	<dir>	output dir
	-ref	<file>	reference file
	-ann	<file>	ann summary file
	-gff	<file>	gff file
	-sv	sv calling default off 
	-cnv	cnv calling default off 
	-RAD	RRL calling default off 
	
	-step	pipeline control
          01 fastq qc
          02 reference prepair
          03 bwa mapping
          04 bam sort
          05 bam mkdup
          06 map-stat
          07 haplotype
          08 gvcf typing
          09 vcf-filter
          10 annovar
          11 sv call
          12 cnv call
          13 variant stat
          14 report
	-stop	pipeline control

  -h         Help

USAGE
        print $usage;
        exit;
}
