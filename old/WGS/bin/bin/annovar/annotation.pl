 #!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fexonic,$fvariant,$outfile,$statfile,$indilist);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"exonic:s"=>\$fexonic,
	"variant:s"=>\$fvariant,
	"indi:s"=>\$indilist,
	"outfile:s"=>\$outfile,
			) or &USAGE;
&USAGE unless ($fexonic and $fvariant and $outfile and $indilist);
my @indi=split(/\,/,$indilist);
open In,$fexonic;
open Out,">$outfile";
print Out "#Chr\tpos\tRef\tAlt\tGpos\tGeneID\tMutType\tMutInfo\tFormat\t",join("\t",@indi),"\n";
my %info;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my (undef,$MutInfo,$MutType,undef,undef,undef,undef,undef,$chr,$pos,undef,$ref,$alt,undef,undef,undef,$GT,@info)=split(/\t/,$_);
	$pos=join("\t",$chr,$pos);
	my @pos;
	foreach my $info (@info) {
		my ($gt,$dep,undef)=split(/\:/,$info);
		push @pos,join(":",$gt,$dep);
	}
	$info{$pos}=join("\t",$MutInfo,$MutType,"GT:AD",@pos);
}
close In;
open In,$fvariant;
my %test;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($postion,$geneID,undef,undef,undef,undef,undef,$chr,$pos,undef,$ref,$alt,undef,undef,undef,$GT,@info)=split(/\t/,$_);
	$pos=join("\t",$chr,$pos);
	if ($postion eq "downstream" || $postion =~ /upstream/ || $postion eq "intergenic" ) {
		$geneID="--";
	}
	$geneID=~s/\(.*\)//g;
	my @pos;
	foreach my $info (@info) {
		my ($gt,$dep,undef)=split(/\:/,$info);
		push @pos,join(":",$gt,$dep);
	}
	if (exists $info{$pos}) {
		print Out join("\t",$pos,$ref,$alt,$postion,$geneID,$info{$pos}),"\n";
	}else{
		print Out join("\t",$pos,$ref,$alt,$postion,$geneID,"--","--","GT:AD",@pos),"\n";
	}
}
close In;
close Out;


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -exonic	<file>	input file name
  -variant	<file>	input file name
  -outfile	<dir>	output dir
  -indi	<str>	indilist,split by ,
  -h         Help

USAGE
        print $usage;
        exit;
}
