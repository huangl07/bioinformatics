#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$split,$protein);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"n:s"=>\$split,
	"p"=>\$protein,
	"dsh"=>\$dsh,
			) or &USAGE;
&USAGE unless ($fIn and $dOut);
my $dNR=" /mnt/ilustre/users/long.huang/DataBase/NR/2017-8-30/nr";
my $dNT="/mnt/ilustre/users/long.huang/DataBase/NT/nt";
my $dKEGG="/mnt/ilustre/users/dna/Environment/biotools/kobas-3.0/seq_pep/ko.pep.fasta ";
my $dGO="/mnt/ilustre/users/long.huang/DataBase/GO/go_weekly-seqdb.fasta";
my $dEggNOG="/mnt/ilustre/users/long.huang/DataBase/EggNOG/new/eggnog4.proteins.all.fa";
my $dSwissprot="/mnt/ilustre/users/long.huang/DataBase/Swissprot/uniref90.fasta";
my $blast="/mnt/ilustre/users/dna/.env//bin/diamond";
if ($protein) {
	$blast.=" blastp";
}else{
	$blast.=" blastn";
}
open In,$fIn;
open SH,">$dsh/alighment.sh";
mkdir $dOut if (!-d $dOut);
$dOut=$ABSOLUTE_DIR($dOut);
$/=">";
$split||=20;
my $n=0;
my %filehand;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	$n++;
	my $hand=$n % $split;
	if (!exists $filehand{$hand}) {
		open $filehand{$hand} ,">$dOut/sub.$n.fasta";
		my $fa="$dOut/sub.$n.fasta";
		my $fname="$dOut/sub.$n";
		print SH "$blast --db $dNR --query $fa --evalue 10e-5 --outfmt 5 --threads 8 --out $fname.nr.blast && ";
		print SH "perl $Bin/bin/NRANNO.pl -i $fname.nr.blast -o $fname.nr.anno && ";
		print SH "$blast --query $fa --db $dSwissprot --out $fname.swissprot.blast --evalue 10e-5 --outfmt 5 --threads 8 &&";
		print SH "perl $Bin/bin/NRANNO.pl -i $fname.swissprot.blast -o $fname.swissprot.anno && ";
		print SH "$blast --query $fa --db $dGO --out $fname.GO.blast --evalue 10e-5 --outfmt 5 --threads 8 && ";
		print SH "perl $Bin/bin/GOANNO.pl -i $fname.GO.blast -o $fname.GO.anno && ";
		print SH "$blast --query $fa --db $dKEGG --out $fname.kegg.blast --evalue 10e-5 --outfmt 6 --threads 8 && ";
		print SH "annotate.py -i $fname.kegg.blast -t blastout:tab -s ko -o $fname.kegg.kobas &&\n";
		print SH "perl $Bin/bin/KEGGanno.pl -i $fname.kegg.kobas -o $fname.kegg.anno"
		print SH "emapper.py -i $dOut/fasta/sub.$n.fasta -m diamond --dmnd_db $dEggNOG  --data_dir $dEggNOG -o $dOut/eggnog/sub.$n.eggnog --cpu 8 \n ";
	}
	print {$filehand{$hand}} ">$_";
}
close In;
close SH;

print SH "perl $Bin/bin/merge_result.pl -i $dOut/result/ -o $dOut/$Key.summary.result -g $fIn";

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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
