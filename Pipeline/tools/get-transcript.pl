#!/mnt/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$gene,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "ref:s"=>\$ref,
        "gff:s"=>\$gff,
		"out:s"=>\$out,
		"gene:s"=>\$gene,
                        ) or &USAGE;
&USAGE unless ($ref and $gff and $gene);
my %stat;
open Ref,$ref;
if($ref=~/gz$/){
        close Ref;
        open Ref,"gunzip -c $ref|";
}
$/=">";
while(<Ref>){
        chomp;
        next if ($_ eq ""|| /^$/);
        my($chr,@seq)=split(/\s+/,$_);
        my $seq=join("",@seq);
        $stat{$chr}{seq}=$seq;
}       
close Ref;

open Gff,$gff;
if($gff=~/gz$/){
	close Gff;
	open Gff,"gunzip -c $gff|";
}
$/="\n";
while (<Gff>){
	chomp;
	next if ($_ eq "" || /^$/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
	my($chr,$source,$type,$start,$end,$score,$strand,$phase,$info,@unkon)=split(/\t/,$_);
	my ($geneid,$transid,$tran,$exonid);
	
	if($type=~/mRNA/){
		#print join("\t",$chr,$type,$start,$end,$strand,$info),"\n";
		$geneid=$1 if($info=~/Parent=(\w+)\;(.*)/);
		$transid=$1 if($info=~/ID=(\w+)\;(.*)/);
		if($geneid eq $gene){
			$tran=$transid;
		}
		next;
	}
	if($type eq /exon/){
		$exonid=$1 if($info=~/ID=(\w+)\;(.*)/);
		$transid=$1 if($info=~/Parent=(\w+)\;(.*)/);
		if($tran eq $transid){
			my $pos=join("\t",$geneid,$transid,$exonid,$start,$end);
			$stat{$chr}{$strand}{type}=$pos;
		}
	}
}
close Gff;

open Out,">$out/$gene.exon.fa";
foreach my $chr(keys %stat ){
	foreach my $strand (keys %{$stat{$chr}}){
		my $seq=$stat{$chr}{seq};
		my($geneid,$trans,$exon,$start,$end)=split(/\t/,$stat{$chr}{$strand}{type});
		my $length= $end - $start + 1;
		my $exonseq=substr($stat{$chr}{seq},$start - 1,$length);
		#print Out ">",$exon,"\n",$exonseq,"\n";
		if($strand=~/\+/){
			print Out ">",$trans,"\:",$exon,"\n",$exonseq,"\n";	
		}else{
			$exonseq=~s/a|A/T/g;
			$exonseq=~s/t|T/A/g;
			$exonseq=~s/g|G/C/g;
			$exonseq=~s/c|C/G/g;
			my @a;
			for(my $i = 0; $i < length($exonseq); $i++){ 
				$a[$i] = substr($exonseq,$i,1);
			}
			@a=reverse@a;
			print Out ">",$exon,"\n";
			print Out join("",@a),"\n";
		}
	}
}
close Out;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : "d,time()-$BEGIN_TIME,"s\n";
########################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        get gene's seq
        eg:
        perl $Script -ref -gff -vcf -out -gene -list

Usage:
  Options:
	-ref	<file>  input ref.fa
	-gff	<file>  input ref.gff
	-out	<file>	output file name
	-gene	<str>	geneid
	-h			Help

USAGE
        print $usage;
        exit;
}
	

