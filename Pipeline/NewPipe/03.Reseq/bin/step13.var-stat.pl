#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$vcflist,$bamlist,$annolist,$svlist,$cnvlist,$dOut,$dShell,$metric,$gff,$chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcflist,
	"anno:s"=>\$annolist,
	"proc:s"=>\$proc,
	"sv:s"=>\$svlist,
	"cnv:s"=>\$cnvlist,
	"gff:s"=>\$gff,
	"chr:s"=>\$chr,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($vcflist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$proc||=20;
$vcflist=ABSOLUTE_DIR("$vcflist");
$annolist=ABSOLUTE_DIR("$annolist");
$dOut=ABSOLUTE_DIR("$dOut");
$dShell=ABSOLUTE_DIR("$dShell");
$svlist=ABSOLUTE_DIR("$svlist") if ($svlist);
$cnvlist=ABSOLUTE_DIR("$cnvlist") if ($cnvlist);
open SH,">$dShell/13.variant-stat.sh";
open In,$vcflist;
my %pop;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($type,$vcf,$annostat)=split(/\s+/,$_);
	if ($type eq "snp") {
		$pop{snp}=$vcf;
		print SH "perl $Bin/bin/snp.stat.pl -i $vcf -o $dOut/snp.stat -m $dOut/snp.matrix && ";
		print SH "Rscript $Bin/bin/diff_matrix.R --i $dOut/snp.matrix --o $dOut/snp.diff\n";
		print SH "perl $Bin/bin/variant_qual.pl -i $vcf -o1 $dOut/snp.depth -o2 $dOut/snp.GQ &&";
		print SH "Rscript $Bin/bin/variant_qual.R --GQ $dOut/snp.GQ --dep $dOut/snp.depth --o $dOut/snp.qual --str SNP\n";
		print SH "perl $Bin/bin/snpEff-stat.pl -i $annostat -o $dOut/snp\n";
	}elsif ($type eq "indel") {
		$pop{indel}=$vcf;
		print SH "perl $Bin/bin/indel.stat.pl -i $vcf -o $dOut/indel.stat -m $dOut/indel.matrix -l $dOut/indel.len && ";
		print SH "Rscript $Bin/bin/indel_len.R --i $dOut/indel.len --o $dOut/ && ";
		print SH "Rscript $Bin/bin/diff_matrix.R --i $dOut/indel.matrix --o $dOut/indel.diff\n";
		print SH "perl $Bin/bin/variant_qual.pl -i $vcf -o1 $dOut/indel.depth -o2 $dOut/indel.GQ && ";
		print SH "Rscript $Bin/bin/variant_qual.R --GQ $dOut/indel.GQ --dep $dOut/indel.depth --o $dOut/indel.qual --str INDEL\n";
		print SH "perl $Bin/bin/snpEff-stat.pl -i $annostat -o $dOut/snp\n";
	}
}
close In;
if ($svlist) {
	open In,$svlist;
	my $max=0;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		my ($sampleID,$sv)=split(/\s+/,$_);
		my $line=`wc -l $sv`;
		chomp $line;
		if ($max < $line) {
			$max=$line;
			$pop{sv}=$sv;
		}
		print SH "perl $Bin/bin/sv.stat.pl -i $sv -o $dOut/$sampleID && ";
		print SH "perl $Bin/bin/varintlen.R -i $dOut/$sample.sv.len -o $dOut/$sample.sv.len\n";
	}
	close In;
}
if ($cnvlist) {
	open In,$cnvlist;
	my $max=0;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		my ($sampleID,$cnv)=split(/\s+/,$_);
		my $line=`wc -l $cnv`;
		chomp $line;
		if ($max < $line) {
			$max=$line;
			$pop{cnv}=$cnv;
		}
		print SH "perl $Bin/bin/cnv.stat.pl -i $cnv -o $dOut/$sampleID && ";
		print SH "perl $Bin/bin/varintlen.R -i $dOut/$sample.cnv.len -o $dOut/$sample.cnv.len\n";
	}
	close In;
}
close In;
print SH "perl $Bin/bin/draw.circos.pl --windows 100000 --snp $pop{snp} --indel $pop{indel} --chrlist $chr --gff $gff --outdir $dOut";
print SH "--sv $pop{sv} " if($svlist);
print SH "--cnv $pop{cnv} " if($cnvlist);
print SH "\n";

close SH;

my $job="perl /mnt/ilustre/users/dna/.env//bin/qsub-sge.pl  --Resource mem=20G --CPU 1 --maxjob $proc  $dShell/13.variant-stat.sh";
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
  -vcf	<file>	input vcf list
  -sv	<file>	input sv list
  -cnv	<file>	input cnv list
  -anno	<file>	input anno list
  -chr	<file>	input chr list
  -gff	<file>	input gff file
  -out	<dir>	output dir
  -proc <num>   number of process for qsub,default 20
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
