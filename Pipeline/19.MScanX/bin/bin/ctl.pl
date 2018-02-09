#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fchr1,$fchr2,$dOut,$fdsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"c1|chrlist1:s"=>\$fchr1,
	"c2|chrlist2:s"=>\$fchr2,
	"o|out:s"=>\$dOut,
			); 		#;/or &USAGE;
&USAGE unless ($fchr1 and $fchr2 and $dOut);
open In1,$fchr1;
open In2,$fchr2;
mkdir $dOut if (!-d $dOut);
mkdir $fdsh if (!-d $fdsh);
my $dir=&DIR($dOut);
my (@chr1,@chr2);
while (<In1>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($before,$after)=split /\s+/,$_;
	push @chr1,$after;
}
close In1;
while(<In2>){
	chomp;
	next if ($_ eq "" || /^$/);
	my ($before,$after)=split /\s+/,$_;
	push @chr2,$after;
}
close In2;
if((scalar @chr1 || scalar @chr2)==0){die "please check two column's chrlist file\n";}
open Circle,">$dir/pro.circle.ctl";
print Circle "1000\t//plot width and height (in pixels)\n";
print Circle join(",",@chr1),",",join(",",@chr2),"\t//chromosomes in the circle\n";
close Circle;
open Dot,">$dir/pro.dot.ctl";
print Dot "800\t//dimension (in pixels) of x axis\n";
print Dot "800\t//dimension (in pixels) of y axis\n";
if(scalar @chr1 <= scalar @chr2){
	print Dot join(",",@chr1),"\t//chromosomes in x axis\n";
	print Dot join(",",@chr2),"\t//chromosomes in y axis\n";
}else{
	print Dot join(",",@chr1),"\t//chromosomes in x axis\n";
	print Dot join(",",@chr2),"\t//chromosomes in y axis\n";	
}
close Dot;
open DualSynteny,">$dir/pro.dualsynteny.ctl";
print DualSynteny "700\t//plot width (in pixels)\n";
print DualSynteny "1000\t//plot height (in pixels)\n";
print DualSynteny join(",",@chr1),"\t//chromosomes in the left column\n";
print DualSynteny join(",",@chr2),"\t//chromosomes in the right column\n";
close DualSynteny;
open Bar,">$dir/pro.bar.ctl";
print Bar "800\t//dimension (in pixels) of x axis\n";
print Bar "800\t//dimension (in pixels) of y axis\n";
print Bar join(",",@chr1),"\t//reference chromosomes\n";
print Bar join(",",@chr2),"\t//target chromosomes\n";
close Bar;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o
Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
