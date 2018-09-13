#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$group,$output);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"gro:s"=>\$group,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($vcf and $group and $output);
open In,$group;
my %sex;
my %info;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$sex,$popid)=split(/\s+/,$_);
	$sex{$sex}++;
	$info{$id}=join("\t",$sex,$popid);
}
close In;
open In,$vcf;
my @Indi;
my %genos;
my @head;
push @head,"IND";
push @head,"SEX";
push @head,"POP";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^##/);
	if (/#/) {
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		push @Indi,@indi;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
		my @alles=split(",",join(",",$ref,$alt));
		next if (scalar @alles >2);
		my @format=split(/\:/,$format);
		push @head,"A";
		my %stat;
		for (my $i=0;$i<@geno;$i++) {
			next if (!exists $info{$Indi[$i]});
			my @info=split(/\:/,$geno[$i]);
			for (my $j=0;$j<@format;$j++) {
				if ($format[$j] eq "GT") {
					my @gt=split(/\//,$info[$j]);
					if ($info[$j] eq "./.") {
						push @{$genos{$Indi[$i]}},9;
					}else{
						$stat{2-$gt[0]-$gt[1]}++;
						push @{$genos{$Indi[$i]}},2-$gt[0]-$gt[1];
					}
				}
			}
		}
		if (scalar keys %stat ==1) {
			pop(@head);
			foreach my $indi (keys %genos) {
				my $test=pop(@{$genos{$indi}})
			}
		}
	}
}
close In;
open Out,">$output";
$sex{F}||=0;
$sex{M}||=0;
my $F=$sex{F};
my $M=$sex{M};
my $sexratio=1;
if ($sex{F} !=0) {
	$sexratio=$M/$F;
}
print Out "insert title here <NM=$sexratio"."NF>\n";
print Out join("\t",@head),"\n";
foreach my $id (sort keys %info) {
	next if (!exists $genos{$id});
	print Out join("\t",$id,$info{$id},@{$genos{$id}}),"\n";
}
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
	"vcf:s"=>\$vcf,
	"gro:s"=>\$group,
	"output:s"=>\$output,
  -h         Help

USAGE
        print $usage;
        exit;
}
