#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dIn,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$dIn,
	"o:s"=>\$out,
			) or &USAGE;
&USAGE unless ($dIn and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input qtl's trit gene's total file's dir 
  -o	<file>	output qtl's gene result xls
  -h         Help

USAGE
        print $usage;
        exit;
}
my @gene=glob("$dIn/*.gene.total");
open Out,">$out";
print Out"#Trait\tChr\tStart\tEnd\tGene Number\tGene with Eff Variation\n";
my %stat;
for (@gene){
	my $trit=(split(/\//,$_))[-1];
	$trit=~s/\.gene\.total//g;
	open In,$_;
	while (<In>){
		chomp;
		if(/^@/){ 
			my($chr,$start,$end,$total,$eff,@undi)=split(/\t/,$_);
			$chr=(split(/\@chr/,$chr))[-1];
			$stat{$trit}{$chr}{anno}=join("\t",$chr,$start,$end,$total,$eff);
		}
	}
	close In;
}
foreach my $trit (keys %stat){
	foreach my $nr (sort{$a<=>$b}keys %{$stat{$trit}}){
		print Out "$trit\t$stat{$trit}{$nr}{anno}\n";
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR
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

