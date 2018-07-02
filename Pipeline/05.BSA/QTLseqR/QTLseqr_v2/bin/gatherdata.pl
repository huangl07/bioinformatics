#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut/3-14.xls";
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || !/\@chr/);
	my @info=split(/\s+/,$_);
	print Out join("\t",@info[0..4]),"\n";
	if (!/#/){
		$stat{NR}{total}+=$info[5];
		$stat{Uni}{total}+=$info[6];
		$stat{KEGG}{total}+=$info[7];
		$stat{GO}{total}+=$info[8];
		$stat{EGGNOG}{total}+=$info[9];
		$stat{NR}{eff}+=$info[10];
		$stat{Uni}{eff}+=$info[11];
		$stat{KEGG}{eff}+=$info[12];
		$stat{GO}{eff}+=$info[13];
		$stat{EGGNOG}{eff}+=$info[14];
	}
}
close In;
close Out;
open Out,">$fOut/3-15.xls";
print Out "#Database\tTotal\tEFF\n";
print Out join("\t","NR",$stat{NR}{total},$stat{NR}{eff}),"\n";
print Out join("\t","Uni",$stat{Uni}{total},$stat{Uni}{eff}),"\n";
print Out join("\t","KEGG",$stat{KEGG}{total},$stat{KEGG}{eff}),"\n";
print Out join("\t","GO",$stat{GO}{total},$stat{GO}{eff}),"\n";
print Out join("\t","EGGNOG",$stat{EGGNOG}{total},$stat{EGGNOG}{eff}),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}
