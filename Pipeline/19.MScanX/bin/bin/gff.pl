#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fgff,$fchrlist);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"p|ProFasta:s"=>\$fIn,
	"g|gff:s"=>\$fgff,
	"c|chrlist:s"=>\$fchrlist,
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn and $fgff and $fOut and $fchrlist);
open In,$fIn;
# $fKey=basename($fKey);
# $fKey=substr(basename($fKey),0,2);
my %profa;
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ); 
	my $a=$1;
	my ($names,@seq)=split /\n/,$_;
	my $name=(split /\s+/,$names)[0];
	my $seq=join("",@seq);
	$profa{$name}=1;
#	print Out ">$name\n$seq\n";
}
close In;
# close Out;
$/="\n";
open In2,$fgff || die "$!\n";
my %pos; 	# pos=>(pos1,pos2,pos3,...)
my %chr;
my %chrs;
while (<In2>) {
	chomp;
        next if ($_ eq "" || /^$/ || /^#/);
	my @info=split /\t/,$_;
#	next if ($info[2] =~ /mRNA/ || /exon/ || /region/ || /RNA/i || /gene/);
	if($info[2] =~ /CDS/i){
		my $proID;
		if($info[8]=~ /Genbank:([^,;]*)/){
			$proID=$1 if(exists $profa{$1});
		}elsif($info[8]=~ /protein_id=([^,;]*)/){
			$proID=$1 if(exists $profa{$1});
		}elsif($info[8]=~ /Name=([^,;]*)/){
			$proID=$1 if(exists $profa{$1});
		}else{die "protein fasta not match gff\n";}
		push @{$pos{$proID}},($info[3],$info[4]);
		$chr{$proID}=$info[0];
		$chrs{$info[0]}=1;  		# proID => chrID;
	}
}
close In2;
my %rename;
#if($fchrlist){
	my $chrlist=&DIR($fchrlist);
	open Chr,"$chrlist";
	while(<Chr>){
		chomp;
		next if ($_ eq "" || /^$/);
		my ($before,$after)=split /\s+/,$_;
		if(!exists $chrs{$before}){die "before name not match gff file";}
		$rename{$before}=$after;
	}
#}
# print Dumper \%pos;die;
open Gff,">$fOut" || die "$!\n";
foreach (sort keys %chr){
	my $lit=&Little(@{$pos{$_}});
	my $lar=&Large(@{$pos{$_}});
	$lit||=0;
	$lar||=0;
#	if($chr{$_}=~ /\w*([0-9]*\w*)/){$chr{$_}=$1;}
	if(exists $rename{$chr{$_}}){
		print Gff "$rename{$chr{$_}}\t$_\t$lit\t$lar\n";
	}else{
		print Gff "$chr{$_}\t$_\t$lit\t$lar\n";
	}
}
close Gff;
#######################################################################################
# print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub Little{
        my $little=shift @_;
        foreach (@_){
                $little=$_ if($_<=$little);
        }
        return $little;
}
sub Large{
        my $large=shift @_;
        foreach (@_){
                $large=$_ if($_>=$large);
        }
        return $large;
}
sub DIR{
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;
		$dir=`pwd`;
		chomp $dir;
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
sub USAGE {
        my $usage=<<"USAGE";
Contact:        qingmei.cui\@majorbio.com;
Script:			$Script
Description:
	protein fasta gff thanslate to MCScanX's gff
	eg:
	perl $Script -p pro.fasta -g ref.gff -c chrlist -o output.MCScanX.gff.name
Usage:
  Options:
  -p	<file>	input protein fasta name
  -g	<file>	input ref.gff name
  -c 	<file>	input chrlist(before name \\t after name)
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
