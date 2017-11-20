#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$dsh,$maf,$mis,$dep,$gro,$chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"gro:s"=>\$gro,
	"out:s"=>\$out,
	"dsh:s"=>\$dsh,
			) or &USAGE;
&USAGE unless ($vcf and $out and $dsh  );
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
$vcf=ABSOLUTE_DIR($vcf);
my $id="pop";
open SH,">$dsh/step04.tree-generic1.sh";
print SH "perl $Bin/bin/vcf2tree.pl -i $vcf -o $out/$id \n";
open SH,">$dsh/step04.tree-generic2.sh";
print SH "cd $out/ && raxml-ng --all --msa $out/$id.fa  --prefix $out/$id --bs-trees 1000 2> $out/$id.raxmlHPC.log && ";
print SH "ln -s $out/RAxML_bipartitionsBranchLabels.$id.tree $out/$id.ml.nwk && ";
print SH "Rscript $Bin/bin/tree.R --infile $out/RAxML_bipartitionsBranchLabels.$id.tree  --outfile $out/$id.ml.tree --raxml 1 ";
if ($gro) {
	$gro=ABSOLUTE_DIR($gro);
	print SH "--group $gro\n";
}
print SH "export OMP_NUM_THREADS=16 && FastTreeMP $out/$id.fasta >  $out/$id.nj.tree && ";
print SH "Rscript $Bin/bin/tree.R --infile $out/$id.nj.tree --outfile $out/$id.nj.tree ";
if ($gro) {
	$gro=ABSOLUTE_DIR($gro);
	 print SH  "--group $gro\n"
}
close SH;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.tree-generic1.sh --CPU 16 --Resource mem=10G";
`$job`;
$job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.tree-generic2.sh --CPU 16 --Resource mem=10G";
`$job`;

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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf files
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
