#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$seqlist,$dShell,$dOut,$misalist);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"seq:s"=>\$seqlist,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
			) or &USAGE;
&USAGE unless ($seqlist and $dOut and $dShell);
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$dOut=ABSOLUTE_DIR("$dOut");
$seqlist=ABSOLUTE_DIR("$seqlist");
$ref=ABSOLUTE_DIR("$ref");
$dShell=ABSOLUTE_DIR("$dShell");
mkdir "$dOut/blastd" if (!-d "$dOut/blastd/");
mkdir "$dOut/blastd/ref" if (!-d "$dOut/blastd/ref");

my ($outlist,$newseq,$newdb,$newmisa,$sample);
`makeblastdb -in $ref -dbtype nucl -parse_seqids -out $dOut/blastd/ref/ref`;
open SH,">$dShell/ssrdiff-1.sh";
open In,$seqlist;
while (<In>){
	chomp;
	next if (/^$/ or "");
	my ($name,$seq,$misa)=split(/\s+/);
	print SH "blastn -db $dOut/blastd/ref/ref -query $seq -out $dOut/$name-1.out -outfmt 6 -evalue 1e-10 -num_alignments 5 &&";
	print SH " perl $Bin/bin/blastout.pl -i $dOut/$name-1.out -k $misa -o $dOut/$name -key ref &&";
	print SH " perl $Bin/bin/reblast.pl -i $dOut/$name.misa-2.out -k $seq -o $dOut/$name.nomatch.seq &&";
	mkdir "$dOut/blastd/$name/" if (!-d "$dOut/blastd/$name/");
	print SH " makeblastdb -in $dOut/$name.nomatch.seq -dbtype nucl -parse_seqids -out $dOut/blastd/$name/$name\n";
	$outlist.=" $dOut/$name.misa-1.out";
	$newseq.=" $dOut/$name.nomatch.seq";
	$newdb.=" $dOut/blastd/$name/$name";
	$newmisa.=" $dOut/$name.misa-2.out";
	$sample.=" $name";
}
close In;
close SH;
#print "$outlist\n";
my (undef,$out1,$out2,undef)=split(/\s+/,$outlist);
my (undef,$newseq1,$newseq2,undef)=split(/\s+/,$newseq);
my (undef,$newdb1,$newdb2,undef)=split(/\s+/,$newdb);
my (undef,$newmisa1,$newmisa2,undef)=split(/\s+/,$newmisa);
my (undef,$sample1,$sample2,undef)=split(/\s+/,$sample);
my $refmisa1=$newmisa1;
my $refmisa2=$newmisa2;
$refmisa1=~s/\-2/\-1/;
$refmisa2=~s/\-2/\-1/;
#print "$refmisa1\n$refmisa2\n";

open SH,">$dShell/ssrdiff-2.sh";
print SH "blastn -db $newdb2 -query $newseq1 -out $newseq1.out -outfmt 6 -evalue 1e-10 -num_alignments 5 &&";
print SH " perl $Bin/bin/blastout.pl -i $newseq1.out -k $newmisa1 -o $newseq1 -key seq && ";
print SH " cat $refmisa1 $newseq1.misa-1.out > $newseq1.seq.result && cat $refmisa2 $refmisa2 >$newseq2.ref.result &&";
print SH " perl $Bin/bin/ssr-stat.pl -i $newseq1.seq.result -k $newseq2.ref.result -o $dOut/$sample1.final.result\n";

print SH "blastn -db $newdb1 -query $newseq2 -out $newseq2.out -outfmt 6 -evalue 1e-10 -num_alignments 5 &&";
print SH " perl $Bin/bin/blastout.pl -i $newseq2.out -k $newmisa2 -o $newseq2 -key seq && ";
print SH " cat $refmisa2 $newseq2.misa-1.out > $newseq2.seq.result && cat $refmisa1 $refmisa1 >$newseq1.ref.result &&";
print SH " perl $Bin/bin/ssr-stat.pl -i $newseq2.seq.result -k $newseq1.ref.result -o $dOut/$sample2.final.result\n";
close SH;

my $job="perl \$qsub --Queue dna --Resource mem=20G --CPU 8  --Nodes 1 $dShell/ssrdiff-1.sh";
#print "$job\n";
	`$job`;

my $job1="perl \$qsub --Queue dna --Resource mem=10G --CPU 8  --Nodes 1 $dShell/ssrdiff-2.sh";
#print "$job1\n";
	`$job1`;
=pod
else{
	my $outlist;
	open In,$seqlist;
	my @line=<In>;
	chomp($line[0]);
	chomp($line[1]);
	my ($sam1,$seq1,$misa1)=split(/\t/,$line[0]);
	my ($sam2,$seq2,$misa2)=split(/\t/,$line[1]);
#	print "$sam1,$seq1,$misa1=======\n$sam2,$seq2,$misa2=========\n";
	`makeblastdb -in $seq1 -dbtype nucl -parse_seqids -out $dOut/blastd/ref/ref`;
	open SH,">$dShell/srrdiff-1.sh";
	print SH "blastn -db $dOut/blastd/ref/ref -query $seq2 -out $dOut/$sam2.out -outfmt 6 -evalue 1e-10 -num_alignments 5 &&";
	print SH " perl /mnt/ilustre/users/chongqing.shi/pipeline/11.SSR/bin/blastout.pl -i $dOut/$sam2.out -k $misa2 -o $dOut/$sam2.misa.out\n";
	close SH;
	open SH,">$dShell/srrdiff-2.sh";
	print SH "perl /mnt/ilustre/users/chongqing.shi/pipeline/11.SSR/bin/ssr-stat.pl -i $misa1 -k $dOut/$sam2.misa.out -o $dOut/final.diff.misa";
	my $job="perl \$qsub --Queue dna --Resource mem=20G --CPU 8  --Nodes 1 $dShell/srrdiff-1.sh";
#	print "$job\n";
	`$job`;
	my $job1="perl \$qsub --Queue dna --Resource mem=10G --CPU 8  --Nodes 1 $dShell/srrdiff-2.sh";
#	print "$job1\n";
   `$job1`;
	close In;
}
=cut
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
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -ref	<file>	ref.fa
  -seq	<file>	list (sample\\tsequence\\tmisaout)
  -out	<dir>	outdir
  -dsh	<dir>	work_sh
USAGE
        print $usage;
        exit;
}
