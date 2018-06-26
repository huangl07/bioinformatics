#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bf,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bf:s"=>\$bf,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bf and $out);
mkdir $out if (!-d $out);
$bf=ABSOLUTE_DIR($bf);
$out=ABSOLUTE_DIR($out);

open IN,"$bf";
my $line=1;
while(<IN>){
    $_=~s/[\n\r]//g;
    my (undef,@result)=split/\t/,$_;
    for(my $i=0;$i<=$#result;$i++){
        open OUT,">>$out/$i.env.result";
        if($line==1){
            print OUT "CHR\tPOS\tP\n";
            }else{
                print OUT "1\t$line\t$result[$i]\n";
            }
    
        }
     $line++;
    }
close IN;
close OUT;
my @result=glob("$out/*.env.result");
for(@result){
    my $file=basename($_);
    `Rscript $Bin/bin/bayenv.r --infile $_ --table $out/$file.table --outfile $out/$file.final`;
    }
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR {#
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
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -bf	<file>	input bayenv result file
  -out	<out>	output dir

  -h         Help

USAGE
        print $usage;
        exit;
}
