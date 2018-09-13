#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$ref,$gff,$RAD,$step,$stop,$SV,$CNV,$realign,$anno,$type);
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
	print LOG "bwa mapping && sort\n"; my $time=time();
	print LOG "########################################\n";
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $fqlist=ABSOLUTE_DIR("$outdir/01.fastq-qc/fq.list");
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
	print LOG "bam mkdup\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/03.bwa-mapping/bam.list");
	my $job="perl $Bin/bin/step04.bam-mkdup.pl -bam $bamlist -out $outdir/04.bam-mkdup -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 5) {
	print LOG "########################################\n";
	print LOG "bam stat\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/04.bam-mkdup/bam.list");
	my $chr=ABSOLUTE_DIR("$outdir/02.ref-config/ref.chrlist");
	my $metric=ABSOLUTE_DIR("$outdir/04.bam-mkdup/metric.list");
	my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step05.map-stat.pl -bam $bamlist -met $metric -dict $dict -chr $chr -out $outdir/05.map-stat -dsh $outdir/work_sh -proc 20";
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
	print LOG "haplotype\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/04.bam-mkdup/bam.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	#my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step06.haplotyper.pl -bam $bamlist -ref $ref -out $outdir/06.haplotype -dsh $outdir/work_sh -proc 20";
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
	print LOG "gvcf typing\n"; my $time=time();
	print LOG "########################################\n";
	my $gvcf=ABSOLUTE_DIR("$outdir/06.haplotype/gvcf.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $dict=ABSOLUTE_DIR("$outdir/02.ref-config/ref.dict");
	my $job="perl $Bin/bin/step07.gvcf-typing.pl -gvcf $gvcf -ref $ref -out $outdir/08.gvcf-typing -dsh $outdir/work_sh";
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
	print LOG "vcf-filter\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/07.gvcf-typing/pop.variant.vcf");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step08.vcf-filter.pl -vcf $vcf -ref $ref -out $outdir/08.vcf-filter -dsh $outdir/work_sh -proc 20";
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
	print LOG "annovar\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/08.vcf-filter/vcf.list");
	my $con=ABSOLUTE_DIR("$outdir/02.ref-config/snpEff.config");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step09.annovar.pl -con $con -vcf $vcf -ref $ref -anno $anno -out $outdir/09.annovar -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step=12 if ($step ne $stop);
}
if ($SV){
	print LOG "########################################\n";
	print LOG "sv call\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/04.bam-mkdup/bam.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $con=ABSOLUTE_DIR("$outdir/02.ref-config/snpEff.config");
	my $job="perl $Bin/bin/step10.sv-call.pl -bam $bamlist -ref $ref -con $con -out $outdir/10.sv-call -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
}
if ($CNV){
	print LOG "########################################\n";
	print LOG "cnv call\n"; my $time=time();
	print LOG "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$outdir/04.bam-mkdup/bam.list");
	my $ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $con=ABSOLUTE_DIR("$outdir/02.ref-config/snpEff.config");
	my $bed=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa.bed");
	my $job="perl $Bin/bin/step12.cnv-call.pl -bam $bamlist  -gff $gff -out $outdir/11.cnv-call -dsh $outdir/work_sh -proc 20";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
}
if ($step == 12) {
	print LOG "########################################\n";
	print LOG "variant stat\n"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$outdir/09.annovar/anno.list");
	my $chr=ABSOLUTE_DIR("$outdir/02.ref-config/ref.chrlist");
	my $gff=ABSOLUTE_DIR("$outdir/02.ref-config/ref/genes.gff");
	my $job="perl $Bin/bin/step12.var-stat.pl -vcf $vcf -chr $chr -gff $gff -out $outdir/12.variant-stat -dsh $outdir/work_sh -proc 20 ";
	#if($SV){
	#	$job.=" -sv $outdir/10.sv-call/sv.anno.primary.vcf ";
	#}
	#if ($CNV) {
	#	$job.=" -cnv $outdir/11.cnv-call/pop.cnv ";
	#}
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
	print LOG "report generic\n"; my $time=time();
	print LOG "########################################\n";
	my $fqdir=ABSOLUTE_DIR("$outdir/01.fastq-qc");
	my $mapstat=ABSOLUTE_DIR("$outdir/05.map-stat");
	my $varstat=ABSOLUTE_DIR("$outdir/12.variant-stat");
	my $vcf=ABSOLUTE_DIR("$outdir/09.annovar/anno.list");
	my $anno=ABSOLUTE_DIR("$outdir/09.annovar");
	my $job="perl $Bin/bin/step13.report.pl -fastqc $fqdir -mapstat $mapstat -variant $varstat -vcf $vcf -out $outdir/13.report -annotate $anno";
	if ($SV) {
		my $sv=ABSOLUTE_DIR("$outdir/10.sv-call/");
		$job.=" -sv $sv ";
	}
	if ($CNV) {
		my $cnv=ABSOLUTE_DIR("$outdir/11.cnv-call/");
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

	-step	pipeline control
          01 fastq qc
          02 reference prepair
          03 bwa mapping and sort
          04 bam mkdup
          05 map-stat
          06 haplotype
          07 gvcf typing
          08 vcf-filter
          09 annovar
          10 sv call
          11 cnv call
          12 variant stat
          13 report
	-stop	pipeline control

  -h         Help

USAGE
        print $usage;
        exit;
}
