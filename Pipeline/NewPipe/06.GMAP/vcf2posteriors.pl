#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($VCF,$fOut,$Pid,$Mid);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$VCF,
	"o:s"=>\$fOut,
	"PID:s"=>\$Pid,
	"MID:s"=>\$Mid,
			) or &USAGE;
&USAGE unless ($VCF and $fOut and $Pid and $Mid);
open In,$VCF;
my (@out1,@out2,@out3,@out4,@out5,@out6);
open Out,">$fOut";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^##/);
	if (/^#/) {
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\s+/,$_);
		push @out1,join("\t","CHR","POS");
		push @out2,join("\t","CHR","POS");
		push @out3,join("\t","CHR","POS");
		push @out4,join("\t","CHR","POS");
		push @out5,join("\t","CHR","POS");
		push @out6,join("\t","CHR","POS");
		push @out1,"Family";
		push @out2,"$Pid";
		push @out3,"0";
		push @out4,"0";
		push @out5,"1";
		push @out6,"0";
		push @out1,"Family";
		push @out2,"$Mid";
		push @out3,"0";
		push @out4,"0";
		push @out5,"2";
		push @out6,"0";
		for (my $i=0;$i<@indi;$i++) {
			if ($Pid eq $indi[$i] || $Mid eq $indi[$i]) {
				next;
			}else{
				push @out1,"Family";
				push @out2,$indi[$i];
				push @out3,"$Pid";
				push @out4,"$Mid";
				push @out5,"0";
				push @out6,"0";
			}
		}
	}else{
		last;
	}
}
close In;
print Out join("\t",@out1),"\n";
print Out join("\t",@out2),"\n";
print Out join("\t",@out3),"\n";
print Out join("\t",@out4),"\n";
print Out join("\t",@out5),"\n";
print Out join("\t",@out6),"\n";
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
