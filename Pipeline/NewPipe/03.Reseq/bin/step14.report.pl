#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fastqc,$mapstat,$variant,$fvcflist,$fsvlist,$fcnvlist,$fannolist,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"fastqc:s"=>\$fastqc,
	"mapstat:s"=>\$mapstat,
	"variant:s"=>\$variant,
	"vcf:s"=>\$fvcflist,
	"sv:s"=>\$fsvlist,
	"cnv:s"=>\$fcnvlist,
	"out:s"=>\$dOut,
#	"proc:s"=>\$proc,
			) or &USAGE;
&USAGE unless ($mapstat and $variant and $fvcflist  and $dOut);
#$proc||=20;
$mapstat=ABSOLUTE_DIR($mapstat);
$variant=ABSOLUTE_DIR($variant);
mkdir $dOut if (!-d $dOut);
mkdir "$dOut/Figure" if (!-d "$dOut/Figure");
mkdir "$dOut/Table" if (!-d "$dOut/Table");
mkdir "$dOut/Result" if (!-d "$dOut/Result");
my $dData="$dOut/Result";
print "Data package\n";
mkdir "$dData/SNP" if (!-d "$dData/SNP");
mkdir "$dData/INDEL"if (!-d "$dData/INDEL");
open In,$fvcflist;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($type,$vcffile)=split(/\s+/,$_);
	if ($type eq "snp") {
		`ln -s $vcffile $dData/SNP`;
	}else{
		`ln -s $vcffile $dData/INDEL`;
	}
}
close In;
close In;
if ($fsvlist) {
	mkdir "$dData/SV" if (!-d "$dData/SV");
	open In,$fsvlist;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		my ($sample,$sv)=split(/\s+/,$_);
		`ln -s $sv $dData/SV`;
	}
	close In;
}
if ($fcnvlist) {
	mkdir "$dData/CNV" if (!-d "$dData/CNV");
	open In,$fcnvlist;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		my ($sample,$cnv)=split(/\s+/,$_);
		`ln -s $cnv $dData/CNV`;
	}
	close In;
}
my @fqstat=glob("$fastqc/*.stat");
my %fqstat;
foreach my $fqstat (@fqstat) {
	open In,$fqstat;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/!Total/ ||/#/);
		my(undef,$sample,$readnum,$basenum,$a,$t,$g,$c,$n,$GC,$Q30,$Q20,undef)=split(/\s+/,$_);
		my $sampleID=(split(/\:/,$sample))[0];
		$fqstat{$sampleID}{readnum}+=$readnum;
		$fqstat{$sampleID}{basenum}+=$basenum;
		$fqstat{$sampleID}{GC}+=$GC*$basenum;
		$fqstat{$sampleID}{Q30}+=$Q30*$basenum;
	}
	close In;
}

print "Data package Done!";
print "map stat package\n";
my @mapstat=glob("$mapstat/*.result.stat");
my %stat;
foreach my $mapstat (@mapstat) {
	my $sampleID;
	open In,$mapstat;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/);
		if (/^#/) {
			$sampleID=(split(/\t/,$_))[-1];
		}else{
			my ($type,$num)=split(/\t/,$_);
			$stat{$sampleID}{$type}=$num;
		}
	}
	close In;
}
open Out,">$dOut/Table/3-3.xls";
print Out "#sampleID\tclean Reads\tclean base\tGC(%)\tQ30(%)\n";
foreach my $sample (sort keys %fqstat) {
	my @out;
	push @out,$sample;
	push @out,$fqstat{$sample}{readnum};
	push @out,$fqstat{$sample}{basenum};
	push @out,sprintf("%.2f",$fqstat{$sample}{GC}/$fqstat{$sample}{basenum});
	push @out,sprintf("%.2f",$fqstat{$sample}{Q30}/$fqstat{$sample}{basenum});
	print Out join("\t",@out),"\n";
}
close Out;

open Out,">$dOut/Table/3-4.xls";
print Out "#sampleID\tclean Reads\tMapped Reads\tMapped Ratio(%)\tProperly Mapped(%)\tDuplicate Ratio(%)\n";
foreach my $sample (sort keys %stat) {
	my @out;
	push @out,$sample;
	push @out,$stat{$sample}{"total reads"};
	push @out,$stat{$sample}{"mapped reads"};
	push @out,sprintf("%.2f",$stat{$sample}{"mapped reads"}/$stat{$sample}{"total reads"}*100);
	push @out,sprintf("%.2f",$stat{$sample}{"proper reads"}/$stat{$sample}{"total reads"}*100);
	push @out,sprintf("%.2f",$stat{$sample}{"duplicate ratio(%)"});
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$dOut/Table/3-5.xls";
print Out "#sampleID\tcover base\tcoverage(1X)\tcoverage(5X)\taverage depth\n";
foreach my $sample (sort keys %stat) {
	my @out;
	push @out,$sample;
	push @out,$stat{$sample}{"cover base"};
	push @out,$stat{$sample}{"genome coverage(1X)"};
	push @out,$stat{$sample}{"genome coverage(5X)"};
	push @out,$stat{$sample}{"average depth"};
	print Out join("\t",@out),"\n";
}
close Out;
print "map stat package Done\n";

open Out,">$dOut/Table/3-6.xls";
open In,"$variant/snp.stat";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my @out=split(/\t/,$_);
	next if (/Total/);
	print Out join("\t",@out[0..6]),"\n";
}
close Out;
close In;
open Out,">$dOut/Table/3-9.xls";
open In,"$variant/indel.stat";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my @out=split(/\t/,$_);
	print Out join("\t",@out[0..4]),"\n";
}
close Out;
my @svstat=glob("$variant/*.sv.stat");
if (scalar @svstat > 0) {
	my %sv;
	my %sample;
	foreach my $sv (@svstat) {
		my $id=(split(/\.sv\.stat/,basename($sv)))[0];
		$sample{$id}=1;
		open In,$sv;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/ || /^#/);
			my ($type,$total,$gene,$alen)=split;
			$sv{$type}{$id}{total}=$total;
			$sv{$type}{$id}{gene}=$gene;
			$sv{$type}{$id}{alen}=abs(sprintf("%.2f",$alen));
		}
		close In;
	}
	open Out,">$dOut/Table/3-12.xls";
	print Out "#type",join("\t",sort keys %sample),"\n";
		foreach my $type (sort keys %sv) {
			my @out;
			push @out,$type;
			foreach my $id (sort keys %sample) {
				$sv{$type}{$id}{total}||=0;
				$sv{$type}{$id}{gene}||=0;
				$sv{$type}{$id}{alen}||=0;
				push @out,join(",",$sv{$type}{$id}{total},$sv{$type}{$id}{gene},$sv{$type}{$id}{alen});
			}
			print Out join("\t",@out),"\n";
		}
	close Out;
}
my @cnvstat=glob("$variant/*.cnv.stat");
if (scalar @svstat > 0) {
	my %cnv;
	my %sample;
	foreach my $cnv (@cnvstat) {
		my $id=(split(/\.cnv\.stat/,basename($cnv)))[0];
		$sample{$id}=1;
		open In,$cnv;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/ || /^#/);
			my ($type,$total,$gene,$alen)=split;
			$cnv{$type}{$id}{total}=$total;
			$cnv{$type}{$id}{gene}=$gene;
			$cnv{$type}{$id}{alen}=abs(sprintf("%.2f",$alen));
		}
		close In;
	}
	open Out,">$dOut/Table/3-13.xls";
	print Out "#type",join("\t",sort keys %sample),"\n";
		foreach my $type (sort keys %cnv) {
			my @out;
			push @out,$type;
			foreach my $id (sort keys %sample) {
				$cnv{$type}{$id}{total}||=0;
				$cnv{$type}{$id}{gene}||=0;
				$cnv{$type}{$id}{alen}||=0;
				push @out,join(",",$cnv{$type}{$id}{total},$cnv{$type}{$id}{gene},$cnv{$type}{$id}{alen});
			}
			print Out join("\t",@out),"\n";
		}
	close Out;
}
if (-f "$variant/snp.effects") {
	`cp $variant/snp.effects $dOut/Table/3-8.xls`
}
if (-f "$variant/snp.region") {
	`cp $variant/snp.region $dOut/Table/3-7.xls`
}
if (-f "$variant/indel.effects") {
	`cp $variant/indel.effects $dOut/Table/3-10.xls`
}
if (-f "$variant/indel.region") {
	`cp $variant/indel.region $dOut/Table/3-11.xls`
}


mkdir "$dOut/Figure/variant" if (!-d "$dOut/Figure/variant");
mkdir "$dOut/Figure/mapstat" if (!-d "$dOut/Figure/mapstat");
mkdir "$dOut/Figure/fastqc" if (!-d "$dOut/Figure/fastqc");
`ln -s $variant/*.pdf $dOut/Figure/variant`;
`ln -s $variant/*.png $dOut/Figure/variant`;
`ln -s $variant/*.svg $dOut/Figure/variant`;
`ln -s $mapstat/*.pdf $dOut/Figure/mapstat`;
`ln -s $mapstat/*.png $dOut/Figure/mapstat`;
`ln -s $fastqc/fig/*.pdf $dOut/Figure/fastqc`;
`ln -s $fastqc/fig/*.png $dOut/Figure/fastqc`;

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
  -fastqc	<dir>	input fastq name
  -mapstat	<dir>	input mapstat name
  -variant	<dir>	input variant name
  -sv	<file>	sv list
  -cnv	<file>	cnv list
  -vcf	<file>	vcf list
  -out	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
