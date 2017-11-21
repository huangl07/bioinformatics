#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$dsh,$chr,$stop,$step);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$out,
	"chr:s"=>\$chr,
	"step:s"=>\$step,
	"stop:s"=>\$stop
			) or &USAGE;
&USAGE unless ($ref and $gff and $out);
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$ref=ABSOLUTE_DIR($ref);
$chr=ABSOLUTE_DIR($chr);
$gff=ABSOLUTE_DIR($gff);
open Log,">$out/work_sh/anno.$BEGIN_TIME.log";
$step=1;
if ($step == 1) {
	print Log "########################################\n";
	print Log "ref-rename \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.newref.pl -i $ref -g $gff -f $chr -o $out/01.newref";
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
	my $ref=ABSOLUTE_DIR("$out/01.newref/ref.fasta");
	my $job="perl $Bin/bin/step02.enzyme.pl -i $ref -o $out/02.enzyme";
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
	print Log "gene blast\n",my $time=time();
	print Log "########################################\n";
	my $ref=ABSOLUTE_DIR("$out/01.newref/ref.fasta");
	my $gff=ABSOLUTE_DIR("$out/01.newref/ref.gff");
	my $job="perl $Bin/bin/step02.blast.pl -i $ref -gff $gff -o $out/03.blast";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if ($step ==4) {
	print Log "########################################\n";
	print Log "gene blast\n",my $time=time();
	print Log "########################################\n";
	my $ref=ABSOLUTE_DIR("$out/01.newref/ref.fasta");
	my $gff=ABSOLUTE_DIR("$out/01.newref/ref.gff");
	my $job="perl $Bin/bin/step02.blast.pl -i $ref -gff $gff -o $out/03.blast";
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
	tree	pca		structure	kinship
	eg:
	perl $Script -i -o -k -c

Usage:

	-vcf	<file>	input vcf files
	-out	<dir>	output dir
	-gro	<file>	input group list
	-maf	<num>	maf default 0.05
	-mis	<num>	mis default 0.3
	-dep	<num>	default 2 
	-step		pipiline control
	-h			Help

USAGE
        print $usage;
        exit;
}
