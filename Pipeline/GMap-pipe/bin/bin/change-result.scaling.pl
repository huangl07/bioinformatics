#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my ($dIn,$dOut,$ratio);
GetOptions(
	"in:s"=>\$dIn,
	"out:s"=>\$dOut,
	"help|?"=>\&USAGE,
	)or &USAGE;
&USAGE unless($dIn and $dOut);
#$ratio||=1000000;
my %stat;

open In,"$dIn/total.map";
open Map,">$dOut/total.map";
open Csv,">$dOut/total.csv";
open Mark,">$dOut/total.marker";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/);
	next if ($_=~/group/);
	my ($id,$pos)=split/\t/;
	my $lg=(split(/\D+/,(split(/\_/,$id))[0]))[-1];
	my $nr=(split(/\_/,$id))[1];
	#$stat{$lg}{$nr}++ ;
	$stat{$lg}{$nr}{pos}=$pos;
	$stat{$lg}{$nr}{id}=$id;
}
close In;
open In,"$dIn/total.csv";
$/="\n";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/);
	if($_=~/Genotype/){
		print Csv "$_\n";
	}else{
		my ($id,$lg,$pos,@ann)=split(/,/,$_);
		my $nr=(split(/\_/,$id))[1];
		#$pos=(split(/\_/,$id))[1];
		#$stat{$lg}{$pos}{id}=$id;
		$stat{$lg}{$nr}{ann}=join(",",@ann);
	}
}
close In;
open In,"$dIn/total.marker";
while (<In>){
	chomp;
	if ($_=~/MarkerID/){
		print Mark "$_\n";
		next ;
	}else{
		my ($id,@undi)=split(/\t/,$_);
		my $lg=(split(/\D+/,(split(/\_/,$id))[0]))[-1];
		#my $pos=(split(/\_/,$id))[1];
		my $nr=(split(/\_/,$id))[1];
		$stat{$lg}{$nr}{info}=$_;
	}
}
my $firstpos;
my $number=1;
foreach my $lg(sort {$a<=>$b} keys %stat){
	print Map "group $lg\n";
	foreach my $nr(sort {$a<=>$b}keys %{$stat{$lg}}){
		if ($number ==1){
			$firstpos=$stat{$lg}{$nr}{pos};
			print Map "$stat{$lg}{$nr}{id}\t0\n";
			print Csv $stat{$lg}{$nr}{id},",",$lg,",",0,",",$stat{$lg}{$nr}{ann},"\n";
			#$stat{$lg}{$pos}{info}=~s/$pos/0/;
			print Mark "$stat{$lg}{$nr}{info}\n";
		}else{
			$stat{$lg}{$nr}{pos}=$stat{$lg}{$nr}{pos} - $firstpos ;
			print Map "$stat{$lg}{$nr}{id}\t",sprintf("%.10f",$stat{$lg}{$nr}{pos}),"\n";
			print Csv $stat{$lg}{$nr}{id},",",$lg,",",sprintf("%.10f",$stat{$lg}{$nr}{pos}),",",$stat{$lg}{$nr}{ann},"\n";
			print Mark $stat{$lg}{$nr}{info},"\n";
		}
		$number++ ;
	}
	$number=1;
}
close In;
close Map;
close Csv;
close Mark;
############################################################
sub USAGE{#
	my $usage=<<"USAGE";

USAGE:
	Options:
	-in	<file>	input file name
	-out	<file>	output file name
	-h	help

USAGE
	print $usage;
	exit;
}
