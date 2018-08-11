#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my ($scaffold,$output,$protein);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"scaffold:s"=>\$scaffold,
	"protein:s"=>\$protein,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($scaffold and $output and $protein);
mkdir $output if (!-d $output);
$scaffold=ABSOLUTE_DIR($scaffold);
$output=ABSOLUTE_DIR($output);
mkdir "$output/ref_config" if (!-d "$output/ref_config");
mkdir "$output/geneblast" if (!-d "$output/geneblast");
mkdir "$output/genewise" if (!-d "$output/genewise");

my $scaffoldname=basename($scaffold);
open SH,">$output/step1.sh";
print SH "ln -s $scaffold $output/ref_config && ";
print SH "cp $Bin/alignscore.txt $output && ";
my $formatdb=`which formatdb`;
chomp $formatdb;
my $blastall=`which blastall`;
chomp $blastall;
print SH "cp $formatdb $output && ";
print SH "cp $blastall $output && ";
print SH "cd $output/ref_config/ && formatdb -i $scaffoldname\n";
close SH;
`perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl $output/step1.sh `;

open SH,">$output/step2.sh";
open In,$protein;
$/=">";
my $n=0;
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my $id=(split(/\s+/,$_))[0];
	open Out,">$output/geneblast/$id.protein.fa";
	print Out ">$_\n";
	close Out;
	print SH "cd $output && genblast_v138_linux_x86_64 -p  genblastg  -q $output/geneblast/$id.protein.fa  -t  $output/ref_config/$scaffoldname  -e  1e-5 -g  T  -f  F  -a 0.5 -d 100000  -r  10  -c  0.5  -s  0  -i  15  -x  20 -n 20 -v 2 -h 0 -j 3  -norepair -gff -cdna -pro -o  $output/geneblast/$id.fa.geneblast \n ";
}
close In;
close SH;
`perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl $output/step2.sh `;
open SH,">$output/step3.sh";
my @DNA=glob("$output/geneblast/*_0");
foreach my $DNA (@DNA) {
	my $id=(split(/\.fa/,basename($DNA)))[0];
	print SH "grep rank:1 $DNA|cut -f 2 -d \"|\"|less -S |sed 's/\:/\\t/g'|sed 's/\\.\\./\\t/g'|seqtk subseq $scaffold - > $output/geneblast/$id.dna.fa && ";
	print SH "genewise  -both  $output/geneblast/$id.protein.fa  $output/geneblast/$id.dna.fa  -alg 333  -silent   -pseudo  -divide \"//\" -trans  -cdna   >  $output/geneblast/$id.genewise.out \n";
}
close SH;
######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
######################################################################################


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
  Options:
	"scaffold:s"=>\$input, #scaffold dir
	"protein:s"=>\$protein,#proteind file
	"output:s"=>\$output,#output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
