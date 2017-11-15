#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Key,$Gff,$Step,$Only,$fChr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"c:s"=>\$fChr,
	"g:s"=>\$Gff,
	"step:s"=>\$Step,
	"only"=>\$Only,
	) or &USAGE;
&USAGE unless ($fIn and $dOut and $Key and $Gff);
my $mkdir=1;
$mkdir=(mkdir $dOut) if (!-d $dOut);
die "Error make dir $dOut" if($mkdir == 0);
$mkdir=(mkdir "$dOut/work_sh") if (!-d "$dOut/work_sh");
die "Error make dir $dOut/work_sh" if($mkdir == 0);
$mkdir=(mkdir "$dOut/NewGenome") if (!-d "$dOut/NewGenome");
die "Error make dir $dOut/NewGenome" if($mkdir == 0);
$mkdir=(mkdir "$dOut/Pre-design") if (!-d "$dOut/Pre-design");
die "Error make dir $dOut/Pre-design" if($mkdir == 0);

die "Error make dir $dOut" if($mkdir == 0);
die "Error input Genome $fIn!\n" if (!-f $fIn ) ;
die "Error input Gff $Gff!\n" if (!-f $Gff);
$fIn = Cwd::abs_path($fIn);
$dOut = Cwd::abs_path($dOut);
$Gff = Cwd::abs_path($Gff) ;
$Step ||= 1;
if (-e $fChr) {
	$fChr=Cwd::abs_path($fChr);
}
#my $qsub="$Bin/qsub-sge.pl";
open Log,">$dOut/$Key.genome.config.log";
if ($Step == 1) {
	print Log "########################################\n";
	print Log "Rename Genome!"; my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/GRename.pl -i $fIn -g $Gff -f $fChr -o $dOut/NewGenome -k $Key && ";
	open SH,">$dOut/work_sh/step1.sh";
	print SH $job,"\n";
	close SH;
	$job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dOut/work_sh/step1.sh --Resource mem=3G --CPU 1 ";
	print Log "$job\n";
	my $return=`$job`;
	die "$job" if ($return !~ /Done/);
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Rename Genome Done and elapsed time : ",time()-$BEGIN_TIME,"s\n";
	print Log "########################################\n";

	$Step++ unless($Only);
}
if ($Step == 2) {
	if (!-f "$dOut/NewGenome/$Key.fasta") {
		die "Error in step1,please check!";
	}
	print Log "########################################\n";
	print Log "Pre-design for GBS or RAD"; 
	print Log "########################################\n";
	my @enzyme=("EcoRI","MseI","PstI","TaqaI");
	open SH,">$dOut/work_sh/step2.sh";
	for (my $i=0;$i<@enzyme;$i++) {
		print SH "perl $Bin/bin/Pre-design.pl -i $dOut/NewGenome/$Key.fasta -o $dOut/Pre-design/$enzyme[$i] -k $Key -e1 $enzyme[$i] && ";
		#perl $Bin/bin/EnzymeRad_size_distribution.pl -I $dOut/Pre-design/$enzyme[$i]/$Key.enzyme.detail -O $dOut/Pre-design/$enzyme[$i]/$Key &&"
		print SH "Rscript $Bin/bin/EnzymeCut_Draw.r -i $dOut/Pre-design/$enzyme[$i]/$Key.enzyme.draw -o $dOut/Pre-design/$enzyme[$i]/$Key\_$enzyme[$i].enzyme.draw\n";
		for (my $j=$i+1;$j<@enzyme;$j++) {
			print SH "perl $Bin/bin/Pre-design.pl -i $dOut/NewGenome/$Key.fasta -o $dOut/Pre-design/$enzyme[$i]\-$enzyme[$j] -k $Key -e1 $enzyme[$i] -e2 $enzyme[$j] && ";
			#perl $Bin/bin/EnzymeGbs_size_distribution.pl -I $dOut/Pre-design/$enzyme[$i]\-$enzyme[$j]/$Key.enzyme.detail -O $dOut/Pre-design/$enzyme[$i]\-$enzyme[$j]/$Key && "
			print SH "Rscript $Bin/bin/EnzymeCut_Draw.r -i $dOut/Pre-design/$enzyme[$i]\-$enzyme[$j]/$Key.enzyme.draw -o $dOut/Pre-design/$enzyme[$i]\-$enzyme[$j]/$Key\_$enzyme[$i]\_$enzyme[$j]\n";
		}
	} 
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  $dOut/work_sh/step2.sh --Resource mem=10G --CPU 1 && cat $dOut/Pre-design/*/$Key.pre-design.stat|sort|uniq > $dOut/Pre-design/$Key.enzyme.stat\n";
	print Log "$job\n";
	my $return=`$job`;
	die "$job" if ($return !~ /Done/);
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Pre-design for GBS or RAD Done and elapsed time : ",time()-$BEGIN_TIME,"s\n";
	print Log "########################################\n";

	$Step++ unless($Only);
}
if ($Step == 3) {
	print Log "########################################\n";
	print Log "Genome Annotation"; 
	print Log "########################################\n";
	my $job="perl $Bin/bin/Annotation/Annotation.pl -i $dOut/NewGenome/$Key.gene.fasta -t nuc -o $dOut/Annotation -k $Key -n 50";
	print Log $job,"\n";
	my $return=`$job`;
	die "$job" if ($return !~ /Done/);
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Genome Annotation Done and elapsed time : ",time()-$BEGIN_TIME,"s\n";
	print Log "########################################\n";
	$Step++ unless($Only);
}




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	Genome configure for all DNA pipeline
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir 

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -c	<file>	input chr list file
  -o	<dir>	output dir,
  -k	<str>	output keys of filename,

  -step	<num>	pipeline control
	  1:	rename
	  2:	pre-design for GBS or RAD
	  3:	genome Annotation
  -only		pipeline control

  -h         Help

USAGE
        print $usage;
        exit;
}
