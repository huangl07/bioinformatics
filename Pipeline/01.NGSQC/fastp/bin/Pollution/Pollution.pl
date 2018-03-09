#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$dOut,$Key,$num,$Step);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
#print $Bin;die;
GetOptions(
	"help|?" =>\&USAGE,
	"fqlist:s"=>\$fqlist,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"n:s"=>\$num,
	"step:s"=>\$Step,
			) or &USAGE;
&USAGE unless ($fqlist and $dOut);
$fqlist=Cwd::abs_path($fqlist);
$num||=10000;
mkdir $dOut if (!-d $dOut);
$dOut=Cwd::abs_path($dOut);
mkdir "$dOut/extract" if (!-d "$dOut/extract");
mkdir "$dOut/blast" if (!-d "$dOut/blast");
mkdir "$dOut/result" if (!-d "$dOut/result");
mkdir "$dOut/worksh" if (!-d "$dOut/worksh");
open Log,">$dOut/$Key.Pollution.log";
$Step ||= 1;
my $blast="/mnt/ilustre/users/dna/.env/bin/blastn";
my $dNT="/mnt/ilustre/users/long.huang/DataBase/NT/nt";
if ($Step == 1) {
	open SH,">$dOut/worksh/step1.sh";
	open In,$fqlist;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		my ($id,$fq1,$fq2)=split(/\s+/,$_);
		print SH "perl $Bin/bin/fastq_extract.pl -1 $fq1 -2 $fq2  -o $dOut/extract/$id.extract.fa -n 10000 \n";
	}
	close In;
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --CPU 1  $dOut/worksh/step1.sh";
	print Log "$job\n";
	my $return=`$job`;die "$job" if ($return !~ /Done/);
	print Log "$job\tdone!\n";
	$Step++;
}
if ($Step == 2) {
	my @fa=glob("$dOut/extract/*.fa");
	open SH,">$dOut/worksh/step2.sh";
	foreach my $fa (@fa) {
		my $fname=basename($fa);
		$fname=~s/\.extract.fa//g;
		print SH "export  LD_LIBRARY_PATH=/mnt/ilustre/app/pub/lib/:/mnt/ilustre/app/pub/lib64/ && ";
		print SH "$blast -query $fa -db $dNT -out $dOut/blast/$fname.nt.blast -evalue 10e-3 -outfmt 5 -num_threads 8 && ";
		print SH "perl $Bin/bin/blast_result.pl -i $dOut/blast/$fname.nt.blast -o $dOut/blast/$fname.nt.blast.best -tophit 1&&";
		print SH "perl $Bin/bin/species_stat_by_blast.pl -i  $dOut/blast/$fname.nt.blast.best -o  $dOut/result/$fname.nt.blast.best.stat &&\n";
	}
	close SH;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  $dOut/worksh/step2.sh --CPU 8";
	print Log "$job\n";
	my $return=`$job`;die "$job" if ($return !~ /Done/);
	print Log "$job\tdone!\n";
	$Step++;
}
if ($Step == 3) {
	my @result=glob("$dOut/result/*.stat");
	open Out,">$dOut/$Key.summary.xls";
	foreach my $result (@result) {
		my $samname=basename($result);
		my $name=(split(/\./,$samname))[0];
		open In,$result;
		my $n=0;
		my $tonum=0;
		print Out $name,"\t";
		my %stat;
		while (<In>) {
			chomp;
			my @info;
			next if ($_ eq "" || /^$/ || /^#/);
			@info=split(/\s+/,$_);
			if ($n <5) {
				$stat{$info[0]}+=$info[2];
			}
			$tonum+=$info[2];
		}
		my @species;
		my @percent;
		my @num;
		foreach my $species (sort {$stat{$b}<=>$stat{$a}} keys %stat) {
			push @species,$species;
			push @percent,sprintf("%.2f",100*$stat{$species}/$tonum),"%";
			push @num,$stat{$species};
			last if (scalar @species == 5);
		}
		print Out join(":",@species),"\t",join(":",@percent),"\t",join(":",@num),"\n";
		close In;
	}
	close Out;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description: check polluction by NT database
	eg:
	perl $Script -i ./ -o  Pollution -k test -n 10000

Usage:
  Options:
  -i	<dir>	input fq dir	gz or not ,both ok
  -o	<dir>	output dir
  -k	<str>	output keys of filename
  
  -n	<num>	fq num by random extract
  
  -step	<num>	pipeline control
		-1	fq extract
		-2	blast 2 nt
		-3	merge result
  -h         Help

USAGE
        print $usage;
        exit;
}
