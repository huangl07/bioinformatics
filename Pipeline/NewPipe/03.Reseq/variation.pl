#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$outdir,$ref,$gff,$RAD,$step,$stop,$SV,$CNV,$realign);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"outdir:s"=>\$outdir,
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
mkdir $outdir if (!-d $outdir);
$ref=ABSOLUTE_DIR($ref);
$gff=ABSOLUTE_DIR($gff);
$fqlist=ABSOLUTE_DIR($fqlist);
$outdir=ABSOLUTE_DIR($outdir);
my ($bamlist,$chrlist,$gvcflist,$vcflist,$svlist,$cnvlist,$dict,$annolist,$metric,$dictfile);
mkdir "$outdir/work_sh" if (!-d "$outdir/work_sh");
$step||=1;
$stop||=-1;
open LOG,">$dOut/work_sh/Resequence.log";
my $prestep=$step;
if ($step == 1) {
	print LOG "########################################\n";
	print LOG "fastq qc\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step01.fastq-qc.pl -fqlist $fqlist -outdir $outdir/01.fastq-qc -dsh $outdir/work_sh";
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
	my $job="perl $Bin/bin/step02.ref-config.pl -ref $ref -gff $gff  -outdir $outdir/02.ref-config -dsh $outdir/work_sh";
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
	$ref=ABSOLUTE_DIR("$outdir/02.ref-config/ref.fa");
	my $job="perl $Bin/bin/step03.bwa-mapping.pl -ref $ref -fqlist $fqlist  -outdir $outdir/03.bwa-mapping -dsh $outdir/work_sh";
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
	$bamlist=ABSOLUTE_DIR("$outdir/03.bwa-mapping/bam.list");
	my $job="perl $Bin/bin/step04.bam-sort.pl -bam $bamlist -outdir $outdir/04.bam-sort -dsh $outdir/work_sh";
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
	print LOG "bam mkdup\n"; my $time=time();
	print LOG "########################################\n";
	$bamlist=ABSOLUTE_DIR("$outdir/04.bwa-sort/bam.list");
	my $job="perl $Bin/bin/step04.bam-mkdup.pl -bam $bamlist -outdir $outdir/05.bam-mkdup -dsh $outdir/work_sh";
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
	print LOG "bam mkdup\n"; my $time=time();
	print LOG "########################################\n";
	$bamlist=ABSOLUTE_DIR("$outdir/05.bwa-mkdup/bam.list");
	$chr=ABSOLUTE_DIR("$outdir/02.ref-config/chr.list");

	my $job="perl $Bin/bin/step06.bam-stat.pl -bam $bamlist -outdir $outdir/05.bam-mkdup -dsh $outdir/work_sh";
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
