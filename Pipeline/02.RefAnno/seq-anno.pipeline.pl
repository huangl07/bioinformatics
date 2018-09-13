#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fa,$out,$dsh,$chr,$stop,$step,$gff,$ref,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"out:s"=>\$out,
	"type:s"=>\$type,
	"step:s"=>\$step,
	"stop:s"=>\$stop
			) or &USAGE;
&USAGE unless ($fa and $out);
mkdir $out if (!-d $out);
$type||="nucl";
$fa=ABSOLUTE_DIR($fa);
$out=ABSOLUTE_DIR($out);
$stop||=-1;
mkdir "$out/work_sh" if (!-d "$out/work_sh");
open Log,">$out/work_sh/anno.$BEGIN_TIME.log";
$step||=2;
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
	$step++ if($stop ne $step);
}
if ($step == 2) {
	print Log "########################################\n";
	print Log "enzyme cut\n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step02.split-fa.pl -fa $fa -out $out/02.split -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step) ;
}
if ($step == 3) {
	print Log "########################################\n";
	print Log "NR ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step03.nr-anno.pl -fa $fa -out $out/03.NR -dsh $out/work_sh -type $type" ;
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
}
if ($step == 4) {
	print Log "########################################\n";
	print Log "KEGG ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step04.kegg-anno.pl -fa $fa -out $out/04.KEGG -dsh $out/work_sh -type $type";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
}
if ($step == 5) {
	print Log "########################################\n";
	print Log "GO ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step05.go-anno.pl -fa $fa -out $out/05.GO -dsh $out/work_sh -type $type";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
}
if ($step == 6) {
	print Log "########################################\n";
	print Log "Uniref ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step06.uniref-anno.pl -fa $fa -out $out/06.Uniref -dsh $out/work_sh -type $type";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
}
if ($step == 7) {
	print Log "########################################\n";
	print Log "EGGNOG ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $job="perl $Bin/bin/step07.eggnog-anno.pl -fa $fa -out $out/07.Eggnog -dsh $out/work_sh -type $type";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
}
if ($step == 8) {
	print Log "########################################\n";
	print Log "merge ANNO\n",my $time=time();
	print Log "########################################\n";
	my $fa=ABSOLUTE_DIR("$out/02.split/fasta.list");
	my $nr=ABSOLUTE_DIR("$out/03.NR/NR.anno");
	my $kegg=ABSOLUTE_DIR("$out/04.KEGG/KEGG.anno");
	my $go=ABSOLUTE_DIR("$out/05.GO/GO.anno");
	my $uniprot=ABSOLUTE_DIR("$out/06.Uniref/Uni.anno");
	my $eggnog=ABSOLUTE_DIR("$out/07.Eggnog/EGGNOG.anno");
	my $job="perl $Bin/bin/step08.merge-result.pl -nr $nr -kegg $kegg -go $go -eggnog $eggnog -uniprot $uniprot  -out $out/08.result\n";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if($stop ne $step);
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
		warn "Warning just for file and dir \n$in\n";
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

	-fa	<file>	input protein fa files
	-out	<dir>	output dir
	-step		pipeline control
	-stop		pipeline control
	-type	sequence type
	-h			Help

USAGE
        print $usage;
        exit;
}
