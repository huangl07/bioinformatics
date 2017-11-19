#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$split);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"n:s"=>\$split,
	"dsh"=>\$dsh,
			) or &USAGE;
&USAGE unless ($fIn and $dOut);
my $dNR=" /mnt/ilustre/users/long.huang/DataBase/NR/2017-8-30/nr";
my $dKEGG="/mnt/ilustre/users/dna/Environment/biotools/kobas-3.0/seq_pep/ko.pep.fasta ";
my $dGO="/mnt/ilustre/users/long.huang/DataBase/GO/go_weekly-seqdb.fasta";
my $dEggNOG="";
my $dSwissprot="/mnt/ilustre/users/long.huang/DataBase/Swissprot/uniref90.fasta";
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
		print SH "diamond blastn --db $dNR --query $fa --evalue 10e-10 --outfmt 5 --threads 8 --out $dOut/sub.$n.fasta.nr.blast\n";
		print SH "diamond blastn --db $dNR --query $fa --evalue 10e-10 --outfmt 5 --threads 8 --out $dOut/sub.$n.fasta.kegg.blast\n";
		print SH "diamond blastn --db $dNR --query $fa --evalue 10e-10 --outfmt 5 --threads 8 --out $dOut/sub.$n.fasta.go.blast\n";
		print SH "diamond blastn --db $dNR --query $fa --evalue 10e-10 --outfmt 5 --threads 8 --out $dOut/sub.$n.fasta.swissprot.blast\n";
		print SH "diamond blastn --db $dNR --query $fa --evalue 10e-10 --outfmt 5 --threads 8 --out $dOut/sub.$n.fasta.eggNog.blast\n";
	}
	print {$filehand{$hand}} ">$_";
}
close In;
close SH;
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
