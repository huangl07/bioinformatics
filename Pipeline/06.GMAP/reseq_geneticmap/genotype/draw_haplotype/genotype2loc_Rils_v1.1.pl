#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$Grade,$All,$add);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"o:s"=>\$fOut,
				"g:s"=>\$Grade,
				"All"=>\$All,
				"add"=>\$add,
				) or &USAGE;
&USAGE unless ((($fIn and $fOut) || $All) and $Grade);
if ($All) {
	$fIn||="./";
	my @files=glob("$fIn/*.genotype");
	foreach my $key1 (@files) {
		open In,$key1;
		my $names=$key1;
		if (length($key1) > 20) {
			$names=substr($key1,-20);
		}
		my $head=<In>;
		my @indi=split(/\s+/,$head);
		my $indiSum=scalar @indi -2;
		my $read;
			my %marker;
		my $markerSum=0;
		while ($read =<In>) {
			chomp $read;
			next if ($read =~ /^$/ || $read eq "" || $read =~ /ID/);
			my ($name,$type,$line)=split(/\s+/,$read,3);
			next if ($type ne "aaxbb");
			$line=~s/aa/a/g;
			$line=~s/bb/b/g;
			$line=~s/ab/h/g;
			$line=~s/--/-/g;
			$line=~s/a-/d/g;
			$line=~s/-b/c/g;
			$marker{$name}=$line;
			$markerSum++;
		}
		close In;
		open (Out,">$key1.loc") or die $!;
		if (defined $add){
			print Out "name = $names\n";
			print Out "popt = Ri$Grade\n" if ($Grade >2);
			print Out "popt = F2\n" if ($Grade ==2);
			print Out "nloc = $markerSum\n";
			$indiSum += 3 ;
			print Out "nind = $indiSum\n\n";
			foreach my $key2 (keys %marker) {
				print Out $key2,"\t","a\tb\th\t",$marker{$key2},"\n";
			}
		}
		else{
			print Out "name = $names\n";
			print Out "popt = Ri$Grade\n" if ($Grade >2);
			print Out "popt = F2\n" if ($Grade ==2);
			print Out "nloc = $markerSum\n";
			print Out "nind = $indiSum\n\n";
			foreach my $key2 (keys %marker) {
				print Out $key2,"\t",$marker{$key2},"\n";
			}
		}
		close (Out) ;
	}
}else{
	open In,$fIn;
	my $head=<In>;
	my @indi=split(/\s+/,$head);
	my $indiSum=scalar @indi -2;
	my $read;
	my %marker;
	my $markerSum=0;
	my @order;
	while ($read =<In>) {
		chomp $read;
		next if ($read =~ /^$/ || $read eq "" );

		my ($name,$type,$line)=split(/\s+/,$read,3);
		next if ($type ne "aaxbb");
		$line=~s/aa/a/g;
		$line=~s/bb/b/g;
		$line=~s/ab/h/g;
		$line=~s/--/-/g;
		$line=~s/a-/d/g;
		$line=~s/-b/c/g;
		$marker{$name}=$line;
		push @order,$name;
		$markerSum++;
	}
	close In;
	open (Out,">$fOut.loc") or die $!;
	if (length($fOut)>20 ) {
		$fOut=substr($fOut,-20);
	}
	if (defined $add){
		print Out "name = $fOut.loc\n";
		print Out "popt = Ri$Grade\n" if ($Grade >2);
		print Out "popt = F2\n" if ($Grade ==2);
		print Out "nloc = $markerSum\n";
		$indiSum += 3 ;
		print Out "nind = $indiSum\n\n";
		foreach my $key (@order) {
			print Out $key,"\t","a\tb\th\t",$marker{$key},"\n";
		}
	}
	else{
		print Out "name = $fOut.loc\n";
		print Out "popt = Ri$Grade\n" if ($Grade >2);
		print Out "popt = F2\n" if ($Grade ==2);
		print Out "nloc = $markerSum\n";
		print Out "nind = $indiSum\n\n";
		foreach my $key (@order) {
			print Out $key,"\t",$marker{$key},"\n";
		}
	}
	close (Out) ;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact: Huang Long <huangl\@biomarker.com.cn> 

Usage:
  Options:
	-i	<file>	input file
	-o	<file>	key of output file name
	-g	<num>	grade of Population ,must Giving,if (g=2) then type fo population if f2
	-All		do all genotype files ;
	-add		add a b h column before samples
	-h         Help

USAGE
	print $usage;
	exit;
}
