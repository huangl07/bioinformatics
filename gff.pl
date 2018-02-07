#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gff,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gff:s"=>\$gff,
	"out:s"=>\$out
			) or &USAGE;
&USAGE unless ($gff and $out);
open In,$gff;
open Out,">$out";
my %genes;
my %mRNAs;
my %out;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($chr,$source,$type,$start,$end,$info,$flag,$scoure,$feature)=split(/\t/,$_);
	my $id;
	if ($feature =~/ID=(.*?);/) {
		$id=$1;
		if ($type eq "gene") {
			$genes{$id}{tag}=1;
			if($feature=~/Name=(.*?);/){
				$genes{$id}{name}=$1;
			}else{
				$genes{$id}{name}="0";
			}
			#print Out join("\t",0,$start,$end,$type,$flag,$scoure,$id,0,$genes{$id}{name}),"\n";	
		}elsif ($type eq "mRNA") {
			$mRNAs{$id}{tag}=1;
			if($feature=~/Parent=(.*?);/i){
				$mRNAs{$id}{parent}=$1;
			}else{
				$mRNAs{$id}{parent}="0";
			}
			if($feature=~/Name=(.*?);/){
				$mRNAs{$id}{name}=$1;
			}else{
				$mRNAs{$id}{name}="0";
			}
			my $g_id=$mRNAs{$id}{parent};
			#print Out join("\t",$id,$start,$end,$type,$flag,$scoure,$id,$mRNAs{$id}{name},$genes{$g_id}{name}),"\n";	
		}elsif ($feature=~/Parent=(.*?);/i) {
			my($t_id,$g_id)=findParent($1);
			push @{$out{$g_id}{$t_id}{$start}},join("\t",$mRNAs{$t_id}{name},$start,$end,$type,$flag,$scoure,$genes{$g_id}{name});	
		}
	}else{
			#push @{$out{0}{0}{$start}},join("\t",0,$start,$end,$type,$flag,$scoure,0,0,0);	
	}
}
close In;
foreach my $id (sort keys %out) {
	next if ($id eq "0");
	foreach my $tid (sort keys %{$out{$id}}) {
		next if ($tid eq "0");
		foreach my $pos (sort {$a<=>$b} keys %{$out{$id}{$tid}}) {
			print Out join("\n",@{$out{$id}{$tid}{$pos}}),"\n";
		}
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub findParent {
	my ($p_id)=@_;
	my ($t_id, $g_id);
	if ($genes{$p_id}{tag}){
		$t_id="0";
		$g_id=$p_id;
		return ($t_id,$g_id);
	}elsif ($mRNAs{$p_id}{tag}){
		$t_id=$p_id;
		$g_id=$mRNAs{$t_id}{parent};
		return ($t_id,$g_id);
	}else{
		return ("0","0");
	}
}

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
	eg:
	perl $Script -gff -out 

Usage:
  Options:
  -gff	<file>	input file name
  -out	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
