#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$out,$dsh,$maf,$mis,$dep,$gro);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
    "group:s"=>\$group,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($vcf and $out and $group );
my %basechange=("AA"=>"A","TT"=>"T","GG"=>"G","CC"=>"C",
				"GT"=>"K","CG"=>"S","AC"=>"M","AG"=>"R",
				"AT"=>"W","CT"=>"Y","NN"=>"N"
);

open In,$group;
my %group;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	my ($id,$gid)=split(/\s+/,$_);
	push @{$group{$gid}},$id;
}
close In;
my %seq;
my @indi;
my $nbase=0;
my $nloci=0;
my @length;
open In,$vcf;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@sample)=split(/\t/,$_);
		next if (/\.\/\./);
		my @alle=split(/,/,join(",",$ref,$alt));
		next if (scalar @alle > 2 || $alt =~ /\*/);
		my @info=split(/\:/,$format);
		$nbase++;
		if ($nbase % 50 == 0) {
			$nloci++;
		}					
		last if ($nloci ==3);
		$length[$nloci]++;
		for (my $i=0;$i<@indi;$i++) {
			my @ginfo=split(/\:/,$sample[$i]);
			for (my $j=0;$j<@info;$j++) {
				if ($info[$j] eq "GT") {
					my @base=split(/\//,$ginfo[$j]);
					my $base;
					if ($base[0] eq ".") {
						$base="NN";
					}else{
						$base=join("",sort($alle[$base[0]],$alle[$base[1]]));
					}
					if (!exists $basechange{$base}) {
						print $base;die;
					}
					my $seqbase=$basechange{$base};
					$seq{$indi[$i]}{$nloci}.=$seqbase;
				}
			}
		}
	}
}
close In;
open Out,">$out";
my $ngroup=scalar keys %group;
print Out "   $ngroup $nloci\n";
my @out;
foreach my $l (@length) {
	push @out,"(s$l)";
}
print Out join(", ",@out),"\n";
foreach my $gid (sort keys %group) {
	my $nind=scalar @{$group{$gid}};
	print Out $nind,"   ","pop",$gid,"\n";
	for (my $i=0;$i<@length;$i++) {
		foreach my $id (@{$group{$gid}}) {
			print Out $id,"      ",$seq{$id}{$i},"\n";
		}
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	filter snp and retrive fa 
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input vcf files
  -group <file> input group list;split by \\t
  -out	<dir>	output dir

  -h         Help

USAGE
        print $usage;
        exit;
}
