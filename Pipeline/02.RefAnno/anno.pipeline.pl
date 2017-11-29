#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$dsh,$chr,$stop,$step,$gff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"chr:s"=>\$chr,
	"out:s"=>\$out,
	"step:s"=>\$step,
	"stop:s"=>\$stop
			) or &USAGE;
&USAGE unless ($ref and $gff and $out);
mkdir $out if (!-d $out);
$out=ABSOLUTE_DIR($out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$ref=ABSOLUTE_DIR($ref);
if ($chr) {
	$chr=ABSOLUTE_DIR($chr);
}
$gff=ABSOLUTE_DIR($gff);
open Log,">$out/work_sh/anno.$BEGIN_TIME.log";
$step||=1;
if ($step == 1) {
	print Log "########################################\n";
	print Log "ref-rename \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.new-ref.pl -ref $ref -gff $gff -out $out/01.newref -dsh $out/work_sh";
	if ($chr) {
		$job .= " -chr $chr ";
	}
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "enzyme cut\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/01.newref/ref.gene.fa");
	my $job="perl $Bin/bin/step02.split-fa.pl -fa $fa -out $out/02.split -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "NR ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step03.nr-anno.pl -fa $fa -out $out/03.NR -dsh $out/work_sh" ;
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "KEGG ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step04.kegg-anno.pl -fa $fa -out $out/04.KEGG -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 5) {
	print Log "########################################\n";
	print Log "GO ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step05.go-anno.pl -fa $fa -out $out/05.GO -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 6) {
	print Log "########################################\n";
	print Log "Uniref ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step06.uniref-anno.pl -fa $fa -out $out/06.Uniref -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step == 7) {
	print Log "########################################\n";
	print Log "EGGNOG ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step07.eggnog-anno.pl -fa $fa -out $out/07.Eggnog -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}


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

Usage:

	-ref	<file>	input ref files
	-out	<dir>	output dir
	-gff	<file>	input gff file
	-chr	<file>	input chr file
	-step		pipeline control
	-stop		pipeline control
	-h			Help

USAGE
        print $usage;
        exit;
}
