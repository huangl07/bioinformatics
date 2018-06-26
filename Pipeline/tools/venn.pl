#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($poplist,$out,$var);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"pop:s"=>\$poplist,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($out);
open IN,$poplist;
my %pop;
my %anno;
my %head;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/) ;
	my ($popid,$pop)=split(/\s+/,$_);
		open In,$pop;
		my $read;
		while ($read=<In>) {
			chomp;
			next if ($read eq "" || $read=~ /^$/);
			my ($chr,$pos,$ref,$alt,@info) = split(/\s+/,$read);
			if ($read=~/^#/) {
				for(my $i=0;$i<@info;$i++){
					$info[$i]=$popid."_".$info[$i];
				}
				my $head=join("\t",@info[0..$#info-1]);
				$head{$popid}=$head;
				next;
			}
			$pop{$chr}{$pos}{$popid}=join("\t",@info[0..$#info-1]);
			$anno{$chr}{$pos}=join("\t",$ref,$alt,$info[-1]);
		}
		close In;
}
close IN;
my %filehand;
open Out,">$out/data.xls";
my %stat;
foreach my $chr (sort keys %pop) {
	foreach my $pos (sort keys %{$pop{$chr}}) {
		my $filename =join("&", sort keys %{$pop{$chr}{$pos}});
		if (!exists $filehand{$filename}) {
			(open $filehand{$filename},">$out/$filename.result")|| die "error $filename";
			my @out;
			foreach my $out (sort keys %head) {
				push @out,$head{$out};
			}
			print {$filehand{$filename}} join("\t","#chr","pos","ref","alt","ann",@out),"\n";
		}
		$stat{$filename}++;
		my @out;
		foreach my $out (sort keys %{$pop{$chr}{$pos}}) {
			push @out,$pop{$chr}{$pos}{$out};
		}
		print {$filehand{$filename}} join("\t",$chr,$pos,$anno{$chr}{$pos},@out),"\n";
	}
}
foreach my $out (sort keys %stat) {
	print Out $out,"\t",$stat{$out},"\n"
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o

Usage:
  Options:
	"pop:s"=>\$poplist, input pop.list
	"out:s"=>\$out,  output dir

USAGE
        print $usage;
        exit;
}
