#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($mapstat,$variant,$fvcflist,$fsvlist,$fcnvlist,$fannolist,$dOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"mapstat:s"=>\$mapstat,
	"variant:s"=>\$variant,
	"vcf:s"=>\$fvcflist,
	"sv:s"=>\$fsvlist,
	"cnv:s"=>\$fcnvlist,
	"anno:s"=>\$fannolist,
	"out:s"=>\$dOut,
#	"proc:s"=>\$proc,
			) or &USAGE;
&USAGE unless ($mapstat and $variant and $fvcflist and $fannolist and $dOut);
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
open In,$fannolist;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($sample,$type,$vcf)=split(/\s+/,$_);
	if ($type eq "snp") {
		`ln -s $vcf $dData/SNP`;
	}else{
		`ln -s $vcf $dData/INDEL`;
	}
}
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
`ln -s $mapstat/*.depth.pdf $dOut/Figure/`;
`ln -s $mapstat/*.depth.png $dOut/Figure/`;
`ln -s $mapstat/*.genome.coverage.pdf $dOut/Figure/`;
`ln -s $mapstat/*.genome.coverage.png $dOut/Figure/`;
`ln -s $mapstat/*.insert.pdf $dOut/Figure/`;
`ln -s $mapstat/*.insert.png $dOut/Figure/`;
print "map stat package Done\n";

my @vcfstat=glob("$variant/*.stat");
foreach my $vcfstat (@vcfstat) {
	if ($vcfstat =~/snp/) {
		open Out,">$dOut/Table/3-6.xls";
		open In,$vcfstat;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/);
			my @out=split(/\t/,$_);
			next if (/Total/);
			print Out join("\t",@out[0..6]),"\n";
		}
		close Out;
		close In;
	}elsif ($vcfstat =~ /indel/) {
		open Out,">$dOut/Table/3-9.xls";
		open In,$vcfstat;
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/);
			my @out=split(/\t/,$_);
			print Out join("\t",@out[0..4]),"\n";
		}
		close Out;
	}
}
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
			next if ($_ eq "" || /^$/);
			my ($type,$total,$gene)=split;
			$sv{$type}{$id}{total}=$total;
			$sv{$type}{$id}{gene}=$total;
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
				push @out,$sv{$type}{$id}{total};
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
			next if ($_ eq "" || /^$/);
			my ($type,$total,$gene)=split;
			$cnv{$type}{$id}{total}=$total;
			$cnv{$type}{$id}{gene}=$total;
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
				push @out,$cnv{$type}{$id}{total};
			}
			print Out join("\t",@out),"\n";
		}
	close Out;
}
if (-f "$variant/snp.position.xls") {
	`cp $variant/snp.position.xls $dOut/3-14.xls`
}
if (-f "$variant/snp.function.xls") {
	`cp $variant/snp.function.xls $dOut/3-15.xls`
}
if (-f "$variant/indel.position.xls") {
	`cp $variant/indel.position.xls $dOut/3-16.xls`
}
if (-f "$variant/indel.function.xls") {
	`cp $variant/indel.function.xls $dOut/3-17.xls`
}



`ln -s $variant/*.pdf $dOut/Figure`;
`ln -s $variant/*.png $dOut/Figure`;
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
  -mapstat	<dir>	input mapstat name
  -variant	<dir>	input variant name
  -sv	<file>	sv list
  -cnv	<file>	cnv list
  -anno	<file>	anno list
  -vcf	<file>	vcf list
  -out	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}
