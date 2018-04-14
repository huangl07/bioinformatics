#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gwas,$trit,$block,$fout1,$fout2);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gwas:s"=>\$gwas,
    "trit:s"=>\$trit,
    "block:s"=>\$block,
	"out1:s"=>\$fout1,
    "out2:s"=>\$fout2,
			) or &USAGE;
&USAGE unless ($gwas and $trit and $block and $fout1 and $fout2);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	retrive gwas result
	eg:
	perl $Script -i -b -o 

Usage:
  Options:
  -gwas	<file>	input  gwas file
  -block    <file>  input blocks file
  -trit    
  -out1	<dir1>
  -out2 <dir2>
  -h         Help

USAGE
        print $usage;
        exit;
}
my %hash;
####count the line number
open IN,$gwas;
my @line=<IN>;
my $count=@line;
close IN;

close IN;
open IN,$block;

while(<IN>){
    $_=~s/[\n\r]//g;
    my (undef,@block)=split/\s+/,$_;
    my ($chr,$max,$mix)=region(@block);
    for(my $i=0;$i<=$#block;$i++){
        $block[$i]=~s/(\D+)(\d+)\_(\d+)/$2\_$3/g;
        $hash{$block[$i]}=join("\t",$chr,$max,$mix);
        }
    }
close IN;
open IN,$gwas;
open OUT1,">$fout1/$trit.bf1.snp";
open OUT2,">$fout2/$trit.bf2.snp";
open RE1,">$fout1/$trit.bf1.region";
open RE2,">$fout2/$trit.bf2.region";
print OUT1 "#Chr\tpos\tpvalue\n";
print OUT2 "#Chr\tpos\tpvalue\n";
print RE1 "#Chr\tstart\tend\n";
print RE2 "#Chr\tstart\tend\n";
my %p_value;
while(<IN>){
    $_=~s/[\n\r]//g;
    next if /^#/;
    my ($chr,$pos,@pvalue)=split/\t/,$_;
    my $key=$chr."_".$pos;
    for(my $i=0;$i<=$#pvalue;$i++){
        next if $pvalue[$i]=~/^$/;
        next if $pvalue[$i]=~/N/;
        if($pvalue[$i]*$count < 0.05){
            print OUT1 "$chr\t$pos\t$pvalue[$i]\n";
            if(exists $hash{$key}){
                print RE1 "$hash{$key}\n";
                }else{
                    my $pos1=$pos-500000;
                    my $pos2=$pos+500000;
                    if($pos1 < 0){$pos1 =0;}
                    print RE1 "$chr\t$pos1\t$pos2\n";
                    }
            }
        if($pvalue[$i] < 0.00005){
            print OUT2 "$chr\t$pos\t$pvalue[$i]\n";
            if(exists $hash{$key}){
                print RE2 "$hash{$key}\n";
                }else{
                    my $pos1=$pos-500000;
                    my $pos2=$pos+500000;
                    if($pos1 < 0){$pos1 =0;}
                    print RE2 "$chr\t$pos1\t$pos2\n";
                    }
            }
        }
}
close OUT1;
close OUT2;
close RE1;
close RE2;

sub region{
    my @regions=@_;
    my @number;
    my @return;
    my $chr;
    for(@regions){
        $_=~/(.*)\_(\d+)/;
        $chr=$1;
        push(@number,$2);
        }
        $chr=~s/Chr//;
    $return[0]= $chr;
    $return[1] = (sort{$a<=>$b} @number)[0];
    $return[2] = (sort{$b<=>$a} @number)[0];
    return @return;
}
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

