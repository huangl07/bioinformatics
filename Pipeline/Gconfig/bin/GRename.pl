#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Key,$Gff,$Step,$Only,$match);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"g:s"=>\$Gff,
	"f:s"=>\$match
	) or &USAGE;
&USAGE unless ($fIn and $dOut and $Gff and $Key);
my @info;
my $mkdir=1;
my %chro;
open In,$match;
while (<In>) {
	chomp;
	next if($_ eq "" ||/^$/);
	my ($id,$for)=split(/\s+/,$_);
	$chro{$id}=$for;
}
close In;
$mkdir=(mkdir $dOut) if (!-d $dOut);
die "Error make dir $dOut" if($mkdir == 0);
die "Error input Genome $fIn!\n" if (!-f $fIn ) ;
die "Error input Gff $Gff!\n" if (!-f $Gff);
open oFa,">$dOut/$Key.fasta";
open In,$fIn;
my %change;
my %sequence;
$/=">";
my %GStat;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($id,@line)=split(/\n/,$_);
	$id=(split(/\s+/,$_))[0];
##	print $id;die;
	my $newname="";
	if (exists($chro{$id})){
		$GStat{Nchrom}++; 
		$newname=$chro{$id};
	}else{
		$GStat{Nscaff}++;
		$newname="sca$GStat{Nscaff}";
	}
	die	"$id\n$newname\n\n error changed! please contact "if ($newname eq "");
	#print $id,"\n",$newname;die;
	$id=(split(/\s+/,$id))[0];
	$change{$id}=$newname;
	my $seq=join("",@line);
	$seq =~ s/(\w+)/\U$1/;
	$sequence{$newname}=$seq;
	$GStat{Nlen}+=length($seq);
	my $nseq=$seq;
	$GStat{GC}+=($nseq=~s/G|C//g);
	print oFa ">$newname\n",$seq,"\n";
}
close oFa;
close In;
open Gff,$Gff;
open oGff,">$dOut/$Key.gff";
open oGene,">$dOut/$Key.gene.fasta";
$/="\n";
my %Gene;
my %outpos;
while (<Gff>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	#CP002684.1      Genbank region  1       30427671        .       +       .       ID=id0;Dbxref=taxon:3702;Name=1;chromosome=1;ecotype=Columbia;gbkey=Src;mol_type=genomic DNA
	#CP002684.1      Genbank gene    3631    5899    .       +       .       ID=gene0;Dbxref=TAIR:AT1G01010;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC domain containing protein 1,T25K16.1,T25K16_1
	my ($chr,@line)=split(/\t/,$_);
	if (!exists $change{$chr}) {
		die "$chr\nerrror rename format! please check or contact long.huang\@majoribio.com";
	}
	my $newchr=$change{$chr};
	print oGff join("\t",$newchr,@line),"\n";
	my $outpos=join("\t",$chr,$line[3],$line[4]);
	next if (exists $outpos{$outpos});
	$outpos{$outpos}=1;
	if ($line[1] eq "gene" || $line[1] eq "mRNA") {
		my $genename="$newchr\_".$line[2]."\_".$line[3]."\_";
		if($line[7] =~ /Genbank:([^\"\'\;]*)/){
			$genename.=$1;
		}elsif($line[7]=~/Name=([^\"\'\;]*)/) {
			$genename.=$1;
		}elsif ( $line[7]=~/gene=([^\"\'\;]*)/) {
			$genename.=$1;
		}elsif ($line[7]=~/geneID=([^\"\'\;]*)/) {
			$genename.=$1;
		}elsif ($line[7] =~ /gene_name=([^\"\'\;]*)/) {
			$genename.=$1;
		}elsif ($line[7] =~/ID=([^\"\'\;]*)/) {
			$genename.="$1";
		}
		next if (exists $Gene{$genename});
		$Gene{$genename}=1;
		my $seq=substr($sequence{$newchr},$line[2]-1,$line[3]-$line[2]+1);
		$GStat{Glen}+=length($seq);
		print oGene ">$genename\n$seq\n";
	}
}
close Gff;
close oGff;
close oGene;
$GStat{Gnum}=scalar keys %Gene;

open oCLog,">$dOut/$Key.changelog";
print oCLog "#Before:修改前\n";
print oCLog "#After:修改后\n";
print oCLog "#Before\tAfter\n";
foreach my $id (sort keys %change) {
	print oCLog "$id\t",$change{$id},"\n";
}
close oCLog;
open Out,">$dOut/$Key.genome.stat";
my $sum=0;
foreach my $id (sort {length($sequence{$b})<=>length($sequence{$a})} keys %sequence) {
	$sum+=length($sequence{$id});
	if ($sum / $GStat{Nlen} < 0.5) {
		$GStat{N50n}++;
		$GStat{N50}=length($sequence{$id});
	}
	if ($sum / $GStat{Nlen} < 0.7) {
		$GStat{N70n}++;
		$GStat{N70}=length($sequence{$id});
	}
	if ($sum / $GStat{Nlen} < 0.9) {
		$GStat{N90n}++;
		$GStat{N90}=length($sequence{$id});
	}
}
$GStat{Nscaff}||=0;
$GStat{Nchrom}||=0;
$GStat{Nlen}||=0;
$GStat{Glen}||=0;
$GStat{Gnum}||=0;
$GStat{GC}||=0;
$GStat{N50}||=0;
$GStat{N50n}||=0;
$GStat{N70}||=0;
$GStat{N70n}||=0;
$GStat{N90}||=0;
$GStat{N90n}||=0;
print Out "#Genome\tNChromosome\tNScaffold\tTotalLenth\tGeneLenth\tGeneNum\tGC(%)\tN50\tN50num\tN70\tN70num\tN90\tN90num\n";
print Out join("\t",$Key,$GStat{Nchrom},$GStat{Nscaff},$GStat{Nlen},$GStat{Glen},$GStat{Gnum},int($GStat{GC}/$GStat{Nlen}*10000)/100,$GStat{N50},$GStat{N50n},$GStat{N70},$GStat{N70n},$GStat{N90},$GStat{N90n}),"\n";
close Out;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	reformat genome,rename scaffold name at genome fa file and gff file
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir 

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -o	<dir>	output dir,
  -k	<str>	output keys of filename,
  -f	<file>	chromosome list;

  -h         Help

USAGE
        print $usage;
        exit;
}
