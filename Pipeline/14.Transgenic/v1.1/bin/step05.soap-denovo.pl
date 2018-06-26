#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$dsh,$out,$inset_lens);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "inset:s"=>\$inset_lens,
	"fqlist:s"=>\$fqlist,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($fqlist and $out and $dsh);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$fqlist=ABSOLUTE_DIR($fqlist);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open In,$fqlist;
my %sample;
$inset_lens||=400;
while (<In>) {
     chomp;
     next if ($_ eq ""|| /^$/);
     my ($id,$fq1,$fq2)=split(/\s+/,$_);
     if(exists $sample{$id}){
         open Out,">>$out/$id.config";
         $sample{$id}++;
         }else{
             open Out,">$out/$id.config";
             $sample{$id}=1;
             }
             print Out "max_rd_len=151\n";
             print Out "[LIB]\n";
             print Out "avg_ins=$inset_lens\n";
             print Out "reverse_seq=0\n";
             print Out "asm_flags=3\n";
             print Out "rank=$sample{$id}\n";
             print Out "q1=$fq1\n";
             print Out "q2=$fq2\n";
             close Out;

    }
close In;
open SH,">$dsh/step04.soap-denovo.sh";
my @config=glob("$out/*config");
open OUT,">$out/soap.list";
foreach my $con  (@config){
    my $id=(split/\/|\./,$con)[-2];
    print SH "SOAPdenovo-127mer all -s $con -K 63 -R -o $out/$id.denovo -p 20 1> $out/$id.ass.log 2>$out/$id.ass.err\n";
    print OUT "$id\t$out/$id.denovo.scafSeq\n";
    }
close SH;
close OUT;
my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.soap-denovo.sh";
`$job`;

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
	perl $Script -psmcfa  -out -dsh -h

Usage:
  Options:
  -fqlist  <file> fq list file  
  -inset <number> insert number;defult 400
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
