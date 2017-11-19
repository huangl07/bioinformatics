#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($i,$cpu,$mem,$step,$ref,$insert,$only,$haplotype,$RAD,$nt,$nct,$parent,$aln,$single,$gff,$stop);
my ($fqlist,$dOut,$key,$bwamethod);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"dOut:s"=>\$dOut,
	"step:s"=>\$step,
	"stop:s"=>\$stop,
) or &USAGE;
&USAGE unless ($fqlist and $dOut);
my %fastq;
$dOut||="./";
my $mkdir=1;
$mkdir=(mkdir "$dOut") if (!-d "$dOut");
$dOut=ABSOLUTE_DIR($dOut);
die "Error make dir $dOut" if ($mkdir == 0);
$mkdir=(mkdir "$dOut/work_sh") if (!-d "$dOut/work_sh");
die "Error make dir $dOut/work_sh" if ($mkdir == 0);
open LOG,">$dOut/work_sh/STACKS.log";
$step||=0;
$cpu=1;
$mem="3G";
$nt=4;
$nct=8;
$stop||=-1;
if ($step == 0) {
	print LOG "########################################\n";
	print LOG "reference prepair\n"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step00.uniform.pl -fqlist $fqlist -out $dOut/00.uniform -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "reference prepair\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}

if ($step==1) {
	print LOG "########################################\n";
	print LOG "ustacks \n"; my $time=time();
	print LOG "########################################\n";
	my $fqlist=ABSOLUTE_DIR("$dOut/00.uniform/fq.uniform.list");
	my $job="perl $Bin/bin/step01.ustacks.pl -fqlist $fqlist -out $dOut/01.ustacks -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "ustacks\n Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step==2){
	print LOG "########################################\n";
	print LOG "cstacks"; my $time=time();
	print LOG "########################################\n";
	my $ulist=ABSOLUTE_DIR("$dOut/01.ustacks/ustacks.list");
	my $job="perl $Bin/bin/step02.cstacks.pl -ulist $ulist -out $dOut/02.cstacks -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "cstacks Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if($stop != $step);
}
if ($step==3){#bwa merge and sort
	print LOG "########################################\n";
	print LOG "cstacks"; my $time=time();
	print LOG "########################################\n";
	my $ulist=ABSOLUTE_DIR("$dOut/01.ustacks/ustacks.list");
	my $clist=ABSOLUTE_DIR("$dOut/02.cstacks/cstacks.list");
	my $job="perl $Bin/bin/step03.sstacks.pl -ulist $ulist -clist $clist -out $dOut/03.sstacks -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "cstacks Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if($stop != $step);
}
if ($step == 4) {
	print LOG "########################################\n";
	print LOG "correct"; my $time=time();
	print LOG "########################################\n";
	my $ulist=ABSOLUTE_DIR("$dOut/01.ustacks/ustacks.list");
	my $clist=ABSOLUTE_DIR("$dOut/02.cstacks/cstacks.list");
	my $slist=ABSOLUTE_DIR("$dOut/03.sstacks/sstacks.list");
	my $job="perl $Bin/bin/step04.correct.pl -ulist $ulist -clist $clist -slist $slist -out $dOut/04.correct -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "correct Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if($stop != $step);
}
if ($step == 5) {
	print LOG "########################################\n";
	print LOG "correct"; my $time=time();
	print LOG "########################################\n";
	my $job="perl $Bin/bin/step05.genotype.pl -in $dOut/04.correct -out $dOut/05.genotype -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "correct Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if($stop != $step);
}
if ($step == 6) {
	print LOG "########################################\n";
	print LOG "cstacks"; my $time=time();
	print LOG "########################################\n";
	my $vcf=ABSOLUTE_DIR("$dOut/05.genotype/batch_1.vcf");
	my $slist=ABSOLUTE_DIR("$dOut/04.correct/sstacks.list");
	my $job="perl $Bin/bin/step05.stacks-stat.pl -vcf $vcf -out $dOut/05.snpstat -slist $slist -dsh $dOut/work_sh";
	print LOG "$job\n";
	`$job`;
	print LOG "$job\tdone!\n";
	print LOG "########################################\n";
	print LOG "cstacks Done and elapsed time : ",time()-$time,"s\n";
	print LOG "########################################\n";
	$step++ if($stop != $step);
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
		warn "Warning just for file and dir $in\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
	my $usage=<<"USAGE";
Contact:long.huang\@majorbio.com;
Script:		$Script
Description: variation
	eg:
	perl $Script -i  data -ref ref.fasta -l sample.list -al all_sample.list

Usage:
	Options:
	-fqlist	<file>	fastq list,each pair fastq.gz in oneline eg
			sample1	fq1	fq2
	-dOut	<dir>	output file dir,forced
	-step	<num>	variant process 
	-stop	<num>	pipeline control
	-h         Help

USAGE
        print $usage;
        exit;
}
 

