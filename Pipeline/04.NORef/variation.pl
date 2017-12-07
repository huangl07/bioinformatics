#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$method,$step,$stop,$sample);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"outdir:s"=>\$outdir,
	"method:s"=>\$method,
	"sample:s"=>\$sample,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($fqlist and $outdir);
$method||="RAD";
mkdir $outdir if (!-d $outdir);
$fqlist=ABSOLUTE_DIR($fqlist);
$outdir=ABSOLUTE_DIR($outdir);
mkdir "$outdir/work_sh" if (!-d "$outdir/work_sh");
$step||=1;
$stop||=-1;
open LOG,">$outdir/work_sh/STACKS.$BEGIN_TIME.log";
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "fastq qc\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastqqc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh -proc 20";
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
	print LOG "fastq uniform\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step02.uniform.pl -fqlist $fqlist -out $outdir/02.uniform -method $method -dsh $outdir/work_sh  -proc 20";
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
	print LOG "ustacks\n"; my $time=time();
	print LOG "########################################\n";
	my $fqlist=ABSOLUTE_DIR("$outdir/02.uniform/fq.list");
	my $job="perl $Bin/bin/step03.ustacks.pl -fqlist $fqlist -out $outdir/03.ustacks -dsh $outdir/work_sh  -proc 20";
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
	print LOG "cstacks\n"; my $time=time();
	print LOG "########################################\n";
	my $ustacks=ABSOLUTE_DIR("$outdir/03.ustacks/ustacks.list");
	my $job="perl $Bin/bin/step04.cstacks.pl -ulist $ustacks -out $outdir/04.cstacks -dsh $outdir/work_sh ";
	$job.="-sample $sample" if($sample);
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
	print LOG "sstacks\n"; my $time=time();
	print LOG "########################################\n";
	my $cstacks=ABSOLUTE_DIR("$outdir/04.cstacks/cstacks.list");
	my $ustacks=ABSOLUTE_DIR("$outdir/03.ustacks/ustacks.list");
	my $job="perl $Bin/bin/step05.sstacks.pl -ulist $ustacks -clist $cstacks -out $outdir/05.sstacks -dsh $outdir/work_sh -proc 20";
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
	print LOG "correct\n"; my $time=time();
	print LOG "########################################\n";
	my $cstacks=ABSOLUTE_DIR("$outdir/04.cstacks/cstacks.list");
	my $ustacks=ABSOLUTE_DIR("$outdir/03.ustacks/ustacks.list");
	my $sstacks=ABSOLUTE_DIR("$outdir/05.sstacks/sstacks.list");
	my $job="perl $Bin/bin/step06.correct.pl -ulist $ustacks -clist $cstacks -slist $sstacks -out $outdir/06.correct -dsh $outdir/work_sh";
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
	print LOG "rcstacks\n"; my $time=time();
	print LOG "########################################\n";
	my $ustacks=ABSOLUTE_DIR("$outdir/06.correct/ustacks.list");
	my $group=ABSOLUTE_DIR("$outdir/04.cstacks/sample.list");
	my $job="perl $Bin/bin/step07.rcstacks.pl -ulist $ustacks -out $outdir/07.rcstacks -dsh $outdir/work_sh -sample $group";
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
	print LOG "rsstacks\n"; my $time=time();
	print LOG "########################################\n";
	my $ustacks=ABSOLUTE_DIR("$outdir/06.correct/ustacks.list");
	my $cstacks=ABSOLUTE_DIR("$outdir/07.rcstacks/cstacks.list");
	my $job="perl $Bin/bin/step08.rsstacks.pl -ulist $ustacks -clist $cstacks -out $outdir/08.rsstacks -dsh $outdir/work_sh -proc 20";
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
	print LOG "genotype\n"; my $time=time();
	print LOG "########################################\n";
	my $ustacks=ABSOLUTE_DIR("$outdir/06.correct/ustacks.list");
	my $cstacks=ABSOLUTE_DIR("$outdir/07.rcstacks/cstacks.list");
	my $sstacks=ABSOLUTE_DIR("$outdir/08.rsstacks/sstacks.list");
	my $job="perl $Bin/bin/step09.genotype.pl -ulist $ustacks -clist $cstacks -slist $sstacks -out $outdir/09.genotype -dsh $outdir/work_sh ";
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
	print LOG "stacks stat\n"; my $time=time();
	print LOG "########################################\n";
	my $stacks=ABSOLUTE_DIR("$outdir/06.correct/ustacks.list");
	my $vcf=ABSOLUTE_DIR("$outdir/09.genotype/batch_1.vcf");
	my $job="perl $Bin/bin/step10.stackstat.pl -ulist $stacks -vcf $vcf -out $outdir/10.stacksstat -dsh $outdir/work_sh ";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step == 11) {
	print LOG "########################################\n";
	print LOG "stacks stat\n"; my $time=time();
	print LOG "########################################\n";
	my $stacks=ABSOLUTE_DIR("$outdir/10.stacksstat/");
	my $vcf=ABSOLUTE_DIR("$outdir/09.genotype/");
	my $fastqc=ABSOLUTE_DIR("$outdir/01.fastq-qc/");
	my $job="perl $Bin/bin/step11.report.pl -statdir $stacks -vcfdir $vcf -fastqc $fastqc -out $outdir/11.report -dsh $outdir/work_sh ";
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
	-sample	<file>	sample list for cstacks [P,M]
	-method	<method>	GBS or RAD
	-step	pipeline control
          01 fastq qc
          02 fastq uniform lenth
          03 ustacks
          04 cstacks
          05 sstacks
          06 correct
          07 rcstacks
          08 rsstacks
          09 genotype
          10 stat
	-stop	pipeline control

  -h         Help

USAGE
        print $usage;
        exit;
}
