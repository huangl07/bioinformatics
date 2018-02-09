#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fgff1,$fIn2,$fgff2,$dOut,$fchr1,$fchr2,$fdsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "p1|profa1:s"=>\$fIn1,
        "g1|refgff1:s"=>\$fgff1,
        "c1|chr1:s"=>\$fchr1,   # input chr.before after 
        "p2:s"=>\$fIn2,
        "g2:s"=>\$fgff2,
        "c2:s"=>\$fchr2,
        "o:s"=>\$dOut,
	"d|dsh:s"=>\$fdsh,
                        );              #;/or &USAGE;
&USAGE unless ($fIn1 and $fgff1 and $fchr1 and $fIn2 and $fgff2 and $fchr2 and $dOut );
$fIn1=&DIR($fIn1);
$fIn2=&DIR($fIn2);
$fgff1=&DIR($fgff1);
$fgff2=&DIR($fgff2);
$fchr1=&DIR($fchr1);
$fchr2=&DIR($fchr2);
mkdir $dOut if(!-d $dOut);
my $dir=&DIR($dOut);
# make result/work_sh
mkdir "$fdsh" if (!-d "$fdsh");
open SH2,">$fdsh/1.3gff.sh";
	print SH2 "perl $Bin/bin/gff.pl -p $fIn1 -g $fgff1 -o $dir/p1.gff -c $fchr1 && touch $dir/gff2.check.OK \n";
	print SH2 "perl $Bin/bin/gff.pl -p $fIn2 -g $fgff2 -o $dir/p2.gff -c $fchr2 && touch $dir/gff1.check.OK \n ";
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $fdsh/1.3gff.sh";
`$job`;
my $gff1=&DIR("$dir/gff1.check.OK");
my $gff2=&DIR("$dir/gff2.check.OK");
`cd $dir && cat p1.gff p2.gff >p.gff`;
`rm $dir/gff1.check.OK && rm $dir/gff2.check.OK`;
close SH2;
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
Script:                 $Script
Description:
        protein fasta gff thanslate to MCScanX's gff
        eg:
        perl $Script -p pro.fasta -g ref.gff -c chrlist -o output.MCScanX.gff.name
Usage:
  Options:
  -p    <file>  input protein fasta name
  -g    <file>  input ref.gff name
  -c    <file>  input chrlist(before name \\t after name)
  -o    <file>  output file name
  -h         Help

USAGE
        print $usage;
        exit;
}
