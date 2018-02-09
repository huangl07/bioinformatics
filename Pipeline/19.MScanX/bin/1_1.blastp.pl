#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$dOut,$fIn2,$fdsh);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"p1:s"=>\$fIn1,
	"p2:s"=>\$fIn2,
	"o|out:s"=>\$dOut,		# main.pl input dir-MCScan
	"d|dsh:s"=>\$fdsh,		# main.pl input $fout/MCScanX :not need dsh in main.
			); 		
&USAGE unless ($fIn1 && $fIn2 && $dOut && $fdsh);
my $step=1;
`mkdir -p $dOut` if (!-d $dOut);
`mkdir -p $fdsh` if (!-d $fdsh);
$fIn1=&DIR($fIn1);
$fIn2=&DIR($fIn2);
my $dir=&DIR($dOut);
if ($step == 1){
	open SH1,">$fdsh/1.1db.sh";
	`ln -s $fIn1 $dir/p1.fa && ln -s $fIn2 $dir/p2.fa`;
	print SH1 "cd $dir && makeblastdb -in p1.fa -dbtype prot -parse_seqids -out p1 -logfile p1.log\n";
	print SH1 "cd $dir && makeblastdb -in p2.fa -dbtype prot -parse_seqids -out p2 -logfile p2.log\n";
	close SH1;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=3G --CPU 1 $fdsh/1.1db.sh";
	`$job`;
        $step++;
}
if ($step == 2){
	open SH2,">$fdsh/1.2blastp.sh";
# blastp self;
        print SH2 "cd $dir && blastp -db p1 -query $fIn1 -out p11.blast -outfmt 6 -num_threads 8 -evalue 1e-10 -num_alignments 5 && touch p11.check.OK\n";
        print SH2 "cd $dir && blastp -db p2 -query $fIn2 -out p22.blast -outfmt 6 -num_threads 8 -evalue 1e-10 -num_alignments 5 && touch p22.check.OK\n";
#blastp inter
        print SH2 "cd $dir && blastp -db p1 -query $fIn2 -out p12.blast -outfmt 6 -num_threads 8 -evalue 1e-10 -num_alignments 5 && touch p12.check.OK\n";
        print SH2 "cd $dir && blastp -db p2 -query $fIn1 -out p21.blast -outfmt 6 -num_threads 8 -evalue 1e-10 -num_alignments 5 && touch p21.check.OK\n";
	close SH2;
	my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl --Resource mem=5G --CPU 8 $fdsh/1.2blastp.sh";
	`$job`;
#	`cd $dir && rm p1.perf p2.perf p1.phr p2.phr p1.pin p2.pin p1.pog p2.pog p1.psd p2.psd p1.psi p2.psi p1.psq p2.psq`;
	my $p11=&DIR("$dir/p11.check.OK");my $p22=&DIR("$dir/p22.check.OK");
	my $p12=&DIR("$dir/p12.check.OK");my $p21=&DIR("$dir/p21.check.OK");
	`cd $dir && cat p11.blast p22.blast p12.blast p21.blast >p.blast`;
}
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
	perl $Script -p1 fa -p2 fa -o <DIR> -d <DIR>
Usage:
  Options:
  -p1	<file>	input fa name
  -p2	<file>	input fa name
  -o	<DIR>	output DIR
  -d	<DIR>	output DIR
  -h         Help

USAGE
        print $usage;
        exit;
}
