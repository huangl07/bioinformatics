#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fdb,$fqu,$qtype,$dtype,$dOut,$Key,$SplitNum,$Resource);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"q:s"=>\$fqu,
	"d:s"=>\$fdb,
	"qtype:s"=>\$qtype,
	"dtype:s"=>\$dtype,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"n:s"=>\$SplitNum,
	"resource:s"=>\$Resource,
			) or &USAGE;
&USAGE unless ($fqu and $fdb and $dOut and $Key);
$qtype||="nucl";
$dtype||="proc";
$SplitNum||=5;
$fqu=Cwd::abs_path($fqu);
$fdb=Cwd::abs_path($fdb);
mkdir $dOut if (!-d $dOut);
$dOut=Cwd::abs_path($dOut);
mkdir "$dOut/worksh" if (!-d "$dOut/worksh");
mkdir "$dOut/blast" if (!-d "$dOut/blast");
mkdir "$dOut/split" if (!-d "$dOut/split");
mkdir "$dOut/result" if (!-d "$dOut/result");
my $makedb="/mnt/ilustre/users/dna/.env//bin/makeblastdb";
my $blastp="/mnt/ilustre/users/dna/.env//bin/blastp";
my $blastn="/mnt/ilustre/users/dna/.env//bin/blastn";
my $blastx="/mnt/ilustre/users/dna/.env//bin/blastx";
my $tblastn="/mnt/ilustre/users/dna/.env//bin/tblastn";
my $blast;
if ($dtype eq "nucl" and $qtype eq "proc") {
	$blast=$blastx;
}elsif ($dtype eq "nucl" and $qtype eq "nucl") {
	$blast=$blastn;
}elsif ($dtype eq "proc" and $qtype eq "nucl"){
	$blast=$tblastn;	
}elsif ($dtype eq "proc" and $qtype eq "proc") {
	$blast=$blastp;
}
my @nhr=glob("$fdb*.*hr");
if (scalar @nhr ==0) {#step1
	my $job="$makedb  -dbtype $dtype -in $fdb  -input_type fasta";
	open SH,">$dOut/worksh/step1.$Key.makedb.sh";
	print SH $job;
	close SH;
	`/mnt/ilustre/users/dna/.env//bin/qsub-sge.pl $dOut/work_sh/step0.$Key.makedb.sh`;
}else{
		open In,$qtype;
		open SH,">$dOut/worksh/step1.$Key.blast.sh";
		my $nseq=0;
		my %handsh;
		while (<In>) {
			chomp;
			next if ($_ eq ""||/^$/);
			$nseq++;
			my $filehand=$nseq % $SplitNum;
			if (!exists $handsh{$filehand}) {
				open $handsh{$filehand},">$dOut/split/$Key.$filehand.fa";
				print SH "$blast -query $dOut/split/$Key.$filehand.fa -dbtype $dtype -db $fdb -out $dOut/blast/$Key.$filehand.blast -evalue 10e-3 -outfmt 5 -num_threads 8 &&";
				print SH "perl $Bin/blast_result.pl -i $dOut/blast/$Key.$filehand.blast -o $dOut/blast/$$Key.$filehand.blast.best\n";
			}
			print {$handsh{$filehand}} ">$_";
		}
		close In;
		close SH;
		foreach my $key (sort keys %handsh) {
			close $handsh{$key};
		}
		my $job="/mnt/ilustre/users/dna/.env//bin/qsub-sge.pl $dOut/worksh/step1.$Key.blast.sh && cat $dOut/blast/*.best > $dOut/result/$Key.blast.all";
		`$job`;
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	blast fasta file to database
	if nseq > 10000, please split it and blast via qsub
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -q	<file>	input query fasta name
  -d	<file>	input database name
  -qtype	<str>	input query type,[nucl|proc],default nucl
  -dtype	<str>	input database type,[nucl|proc],default nucl
  -o	<dir>	output dir name
  -k	<str>	output keys of filename

  -n	<num>	input split number for qsub default 5
  -resource	<str>	input qsub resource default "mem=50G"

  -h         Help

USAGE
        print $usage;
        exit;
}
