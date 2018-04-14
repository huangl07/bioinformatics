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
  -i	<file>	input qtl's trit vcf file's dir 
  -o	<file>	output qtl's vcf result xls
  -h         Help

USAGE
        print $usage;
        exit;
}
my @vcf=glob("$dIn/*.vcf.total");
open Out,">$out";
print Out"#Trait\tChr\tStart\tEnd\tSNP Number\tEff SNP\tInDel Number\tEff InDel\n";
my %stat;
for (@vcf){
	my $trit=(split(/\//,$_))[-1];
	$trit=~s/\.vcf\.total//g;
	open In,$_;
	while (<In>){
		chomp;
		if(/^@/){ 
			#my($chr,$start,$end,$snp,$effsnp,$indel,$effindel)=split(/\t/,$_);
			my $chr=(split(/\@chr/,(split(/\t/,$_))[0]))[1];
			$stat{$trit}{$chr}{anno}=(split(/\@chr/,$_))[1] ;
		}
	}
	close In;
}
foreach my $trit (keys %stat){
	foreach my $chr (sort{$a<=>$b}keys %{$stat{$trit}}){
		print Out "$trit\t$stat{$trit}{$chr}{anno}\n";
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

