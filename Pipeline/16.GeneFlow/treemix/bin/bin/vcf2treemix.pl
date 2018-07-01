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
my %gid;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	my ($id,$gid)=split(/\s+/,$_);
	$group{$id}=$gid;
	$gid{$gid}=1;
}
close In;
open In,$vcf;
open(Out,"| gzip - > $out");
my @indi;
print Out join("\t",sort keys %gid),"\n";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@sample)=split(/\t/,$_);
		my @alle=split(/,/,join(",",$ref,$alt));
		next if (scalar @alle > 2 || $alt =~ /\*/);
		my @info=split(/\:/,$format);
		my %stat;
		for (my $i=0;$i<@indi;$i++) {
			next if (!exists $group{$indi[$i]});
			my @ginfo=split(/\:/,$sample[$i]);
			for (my $j=0;$j<@info;$j++) {
				next if ($info[$j] ne "GT");
				my @base=split(/\//,$ginfo[$j]);
				next if ($base[0] eq ".");
				if (!defined $base[1]) {
					print Dumper @base;die;
				}
				next if ($base[1] eq ".");
				$stat{$group{$indi[$i]}}{$alle[$base[0]]}++;
				$stat{$group{$indi[$i]}}{$alle[$base[1]]}++;
			}
		}
		my @out;
		foreach my $gid (sort keys %gid) {
			$stat{$gid}{$ref}||=0;
			$stat{$gid}{$alt}||=0;
			push @out,join(",",$stat{$gid}{$ref},$stat{$gid}{$alt});
		}
		print Out join("\t",@out),"\n";
	}
}
close Out;
close In;
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
