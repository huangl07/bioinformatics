#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($id,$fqlist,$dsh,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"id:s"=>\$id,
    "fqlist:s"=>\$fqlist,
    "dsh:s"=>\$dsh,
    "out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($id and $fqlist and $out and $dsh);
mkdir $out if (!-d $out);
mkdir $dsh if (!-d $dsh);
$id=ABSOLUTE_DIR($id);
$fqlist=ABSOLUTE_DIR($fqlist);
$out=ABSOLUTE_DIR($out);
$dsh=ABSOLUTE_DIR($dsh);
open IN,"$fqlist";
my %hash;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($id,$fq1,$fq2)=split/\t/,$_;
    $hash{$id}{fq1}=$fq1;
    $hash{$id}{fq2}=$fq2;
    }
close IN;
open SH,">$dsh/step04.retrive.fq.sh";
open OUT,">$out/insert.fq.list";
my @ids=glob("$id/*.id");
for(@ids){
    my $file=basename($_);
    my $id=(split/\./,$file)[0];
    print $id,"\n";
    print SH "perl $Bin/bin/retrive.fq.pl -i $hash{$id}{fq1} -o $out/$file.b1.fq -list $_ && perl $Bin/bin/retrive.fq.pl -i $hash{$id}{fq2} -o $out/$file.b2.fq -list $_\n";
    print OUT "$file\t$out/$file.b1.fq\t$file\t$out/$file.b2.fq\n";
    }

my $job="perl /mnt/ilustre/users/dna/.env//bin//qsub-sge.pl $dsh/step04.retrive.fq.sh";
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
	perl $Script 

Usage:
  Options:
  -id  <dir> id dir
  -fqlist   <file>  fqlist file
  -out	<dir>	output dir
  -dsh	<dir>	output work shell

  -h         Help

USAGE
        print $usage;
        exit;
}
