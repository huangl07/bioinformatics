#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($trit,$dir,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"trit:s"=>\$trit,
    "dir:s"=>\$dir,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($trit and $dir and $out);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	retrive gwas result
	eg:
	perl $Script -trit -dir -out 

Usage:
  Options:
  -trit
  -dir  <dir>   
  -out	<dir>
  -h         Help

USAGE
        print $usage;
        exit;
}
my %hash;
    if(-e "$dir/GAPIT..$trit.GWAS.Results.csv"){
        open F,"$dir/GAPIT..$trit.GWAS.Results.csv";
        while(<F>){
            $_=~s/[\n\r]//g;
            next if /^SNP/;
            my ($snp,$chr,$pos,$pvalue,@others)=split/\,/,$_;
            $hash{$chr}{$pos}{GAPIT}=$pvalue;
        }
        close F;
    }

     if(-e "$dir/MVP.$trit.FarmCPU.csv"){
         open F,"$dir/MVP.$trit.FarmCPU.csv";
          while(<F>){
              $_=~s/[\n\r]//g;
              next if $_=~/Chrom/;
              my ($snp,$chr,$pos,undef,$pvalue)=split/\,/,$_;
              $chr=~s/"(.*)"/$1/g;$chr=~s/\s+//g;
              $pos=~s/"(.*)"/$1/g;$pos=~s/\s+//g;
              $pvalue=~s/"(.*)"/$1/g;$pvalue=~s/\s+//g;
              $hash{$chr}{$pos}{FarmCPU}=$pvalue;
         }
         close F;
         }
      if(-e "$dir/MVP.$trit.GLM.csv"){
          open F,"$dir/MVP.$trit.GLM.csv";
          while(<F>){
              $_=~s/[\n\r]//g;
              next if $_=~/Chrom/;
              my ($snp,$chr,$pos,undef,$pvalue)=split/\,/,$_;
              $chr=~s/"(.*)"/$1/g;$chr=~s/\s+//g;
              $pos=~s/"(.*)"/$1/g;$pos=~s/\s+//g;
              $pvalue=~s/"(.*)"/$1/g;$pvalue=~s/\s+//g;
              $hash{$chr}{$pos}{GLM}=$pvalue;
              }
               close F;
          }
      if(-e "$dir/MVP.$trit.MLM.csv"){
          open F,"$dir/MVP.$trit.MLM.csv";
          while(<F>){
              $_=~s/[\n\r]//g;
              next if $_=~/Chrom/;
              my ($snp,$chr,$pos,undef,$pvalue)=split/\,/,$_;
              $chr=~s/"(.*)"/$1/g;$chr=~s/\s+//g;
              $pos=~s/"(.*)"/$1/g;$pos=~s/\s+//g;
              $pvalue=~s/"(.*)"/$1/g;$pvalue=~s/\s+//g;
              $hash{$chr}{$pos}{MLM}=$pvalue;
              }
              close F;
          }
       #tassel
      #if(-e "$dir/$trit_tassel2.txt"){
       #   
      #    }
open OUT,">$out/$trit.combine.txt";
print OUT "#chr\tpos\tGAPIT\tMVP_FarmCPU\tMVP_GLM\tMVP_MLM\n";
foreach my $keys (sort keys %hash){
    foreach my $pos (sort keys %{$hash{$keys}}){
        print OUT "$keys\t$pos\t$hash{$keys}{$pos}{GAPIT}\t$hash{$keys}{$pos}{FarmCPU}\t$hash{$keys}{$pos}{GLM}\t$hash{$keys}{$pos}{MLM}\n";
    }
}
close OUT;

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

