#!/mnt/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$vcf,$list,$gene,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "ref:s"=>\$ref,
        "gff:s"=>\$gff,
        "vcf:s"=>\$vcf,
		"out:s"=>\$out,
		"list:s"=>\$list,
		"gene:s"=>\$gene,
                        ) or &USAGE;
&USAGE unless ($ref and $gff and $vcf and $list);
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
	next if ($_ eq "" || /^$/|| /^#/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
	my ($chr,$seqty,$type,$start,$end,undef,undef,undef,$info)=split(/\t/,$_);
	my $geneid;
	if($info=~/ID=([^;]*)/){
		$geneid=$1;
	}
	$stat{$chr}{$geneid}{start}=$start;
	$stat{$chr}{$geneid}{end}=$end;
}
close Gff;

my @samples;
open List,$list;
while(<List>){
	chomp;
	next if($_ eq "" || /^$/);
	push @samples,$_;
}
close List;

my($geneinfo,$geneseq);
foreach my $chr(keys %stat ){
	foreach my $geneid(sort keys %{$stat{$chr}}){
		if($gene eq $geneid){
			my $length=$stat{$chr}{$geneid}{end} - $stat{$chr}{$geneid}{start} + 1;
			$geneinfo=join(",",$chr,$stat{$chr}{$geneid}{start},$stat{$chr}{$geneid}{end});
			$geneseq=substr($stat{$chr}{seq},$stat{$chr}{$geneid}{start} - 1,$length);
			#print $geneinfo,"\t",$geneseq,"\n";
			#die;
		}
	}
}

open In,$vcf;
if($vcf=~/gz$/){
	close In;
	open In,"gunzip -c $vcf|";
}
my @sample;
while (<In>){
	chomp;
	next if ($_ =~ "" || /^$/|| /##/);
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@undi)=split(/\t/,$_);#chr1	1072	chr1_1072	A	C	154.84	PASS	AC=2;AF=1.00;AN=2;ANN=C|upstream_gene_variant|MODIFIER|NAC001|gene0|transcript|rna0|protein_coding||c.-2688A>C|||||2559|,C|intergenic_region|MODIFIER|CHR_START-NAC001|CHR_START-gene0|intergenic_region|CHR_START-gene0|||n.1072A>C||||||;DP=6;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=25.81;SOR=2.303;set=variant	GT:AD:DP:GQ:PL	1/1:0,6:6:18:183,18,0
	if($chr eq "#CHROM"){
		@sample=@undi;
		#print join("\t",@sample),"\n";
	}else{
		my ($genechr,$genestart,$geneend)=split(/\,/,$geneinfo);
		next if($chr ne $genechr);	#!=
		if($pos>=$genestart and $pos<= $geneend){
			my @alt=split(/,/,join(",",$ref,$alt));
			my (%ale,%len);
			for (my $i=0;$i<@alt;$i++) {
				$ale{$alt[$i]}=$i;
				$len{length($alt[$i])}=1;
			}
			my $snptype="SNP";
			$snptype="INDEL"  if (scalar keys %len > 1);

			my $snppos=$pos - $genestart;
			my @format=split(/:/,$format);
			for (my $i=0;$i<@sample;$i++) {
				my $sample=$sample[$i];
				my @info=split(/:/,$undi[$i]);
				for (my $j=0;$j<@info;$j++) {
					if ($format[$j] eq "GT") {
						next if ($info[$j] eq "./." || $info[$j] eq "0/0");
						my ($g1,$g2)=split(/\//,$info[$j]);
						#print "$id\t$alt[$g1]\/$alt[$g2]\t$snppos\n";
						if ($g1 eq $g2) {
							$stat{$sample}{$snppos}{type}=$alt[$g1] ;
						}else{
							my $type=length$ref ;
							if($snptype eq "SNP"){
								my $basetype=join("\/",$alt[$g1],$alt[$g2]);
								$stat{$sample}{$snppos}{$type}="R" if($basetype eq "A\/G");
								$stat{$sample}{$snppos}{$type}="M" if($basetype eq "A\/C");
								$stat{$sample}{$snppos}{$type}="W" if($basetype eq "A\/T");
								$stat{$sample}{$snppos}{$type}="Y" if($basetype eq "C\/T");
								$stat{$sample}{$snppos}{$type}="K" if($basetype eq "G\/T");
								$stat{$sample}{$snppos}{$type}="S" if($basetype eq "G\/C");
							}else{
								my $lang;
								if(length$alt[$g1] >= length$alt[$g2]){
									$lang = $alt[$g1];
								}else{
									$lang = $alt[$g2];
								}
								$stat{$sample}{$snppos}{$type}=$lang;
							}
						}
					}
				}
			}
		}
	}
}
close In;
open Out,">$out/$gene.fa";
foreach my $samples(@samples){
	foreach my $sample (keys %stat){
		if($samples eq $sample){
			#print Out ">$samples\n";
			my $seq=$geneseq;
			
			my ($left,$right,$cutpos);
			foreach my $snppos(sort{$a<=>$b}keys %{$stat{$sample}}){
				foreach my $type(keys %{$stat{$sample}{$snppos}}){
					$cutpos=$snppos + $type;
					$left=substr($seq,0,$snppos);
					$right=substr($seq,$cutpos);
					$seq=join("",$left,$stat{$sample}{$snppos}{$type},$right);
				}
			}
			print Out ">$samples\n$seq\n";
		}
	}
}
close Out;
my $job="clustalo -i $out/$gene.fa -o $gene.diff.phy --outfmt=phy ";
`$job`;

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
	-vcf	<file>	pop.final.vcf
	-out	<file>	output file name
	-gene	<str>	geneid
	-list	<file>	input sample.list
	-h			Help

USAGE
        print $usage;
        exit;
}
	

