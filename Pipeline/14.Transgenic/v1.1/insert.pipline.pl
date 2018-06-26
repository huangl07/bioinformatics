#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$insert,$out,$dsh,$len,$fqlist,$step,$stop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"insert:s"=>\$insert,
    "len:s"=>\$len,
    "fqlist:s"=>\$fqlist,
	"out:s"=>\$out,
	"step:s"=>\$step,
    "stop:s"=>\$stop,
			) or &USAGE;
&USAGE unless ($ref and $insert and $out);
mkdir $out if (!-d $out);
mkdir "$out/work_sh" if (!-d "$out/work_sh");
$out=ABSOLUTE_DIR($out);
$ref=ABSOLUTE_DIR($ref);
$insert=ABSOLUTE_DIR($insert);
$fqlist=ABSOLUTE_DIR($fqlist);
$step||=1;
$stop||=-1;
$len||=400;
open Log,">$out/work_sh/pop.$BEGIN_TIME.log";
if ($step == 1) {
	print Log "########################################\n";
	print Log "step01.index \n",my $time=time();
	print Log "########################################\n";
	my $job="perl $Bin/bin/step01.index.pl -ref $ref -insert $insert -out $out/step01.index -dsh $out/work_sh";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "step01.index Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==2) {
	print Log "########################################\n";
	print Log "bwa-mappping \n",my $time=time();
	print Log "########################################\n";
    my $ref=ABSOLUTE_DIR("$out/step01.index/pop.fa");
	my $job="perl $Bin/bin/step02.bwa-mapping.pl  -ref $ref -fqlist $fqlist -out $out/step02.bwa-mapping -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "bwa-mappping Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
    $step++ if ($step ne $stop);
}
if ($step ==3) {
	print Log "########################################\n";
	print Log "bam-insert \n",my $time=time();
	print Log "########################################\n";
	my $bamlist=ABSOLUTE_DIR("$out/step02.bwa-mapping/bam.list");
	my $job="perl $Bin/bin/step03.bam-insert.pl -bamlist $bamlist -out $out/step03.bam-insert -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "bam-insert Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ if ($step ne $stop);
}
if ($step ==4){
    print Log "########################################\n";
    print Log "retrive-fq \n",my $time=time();
    print Log "########################################\n";
    my $job="perl $Bin/bin/step04.retrive.fq.pl -id $out/step03.bam-insert/ -id $out/step03.bam-insert/ -fqlist $fqlist -out $out/step04.retrive-fq -dsh $out/work_sh ";
    `$job`;
    print Log "$job\tdone!\n";
    print Log "########################################\n";
    print Log "bam-insert Done and elapsed time : ",time()-$time,"s\n";
    print Log "########################################\n";
    $step++ if ($step ne $stop);
    }
if ($step ==5) {
	print Log "########################################\n";
	print Log "soap-denovo \n",my $time=time();
	print Log "########################################\n";
    my $fq=ABSOLUTE_DIR("$out/step04.retrive-fq/insert.fq.list");
	my $job="perl $Bin/bin/step05.soap-denovo.pl -fqlist $fq -out $out/step05.soap-denovo/ -dsh $out/work_sh ";
	print Log "$job\n";
	`$job`;
	print Log "$job\tdone!\n";
	print Log "########################################\n";
	print Log "soap-denovo Done and elapsed time : ",time()-$time,"s\n";
	print Log "########################################\n";
	$step++ ;
}
if($step ==6){
    print Log "########################################\n";
    print Log "blast \n",my $time=time();
    print Log "########################################\n";
    my $list=ABSOLUTE_DIR("$out/step05.soap-denovo/soap.list");
    my $db="$out/step01.index/dbname";
    my $job="perl $Bin/bin/step06.blast.pl -list $list -db $db -out $out/step06.blast -dsh $out/work_sh";
    print Log "$job\n";
    `$job`;
    print Log "$job\tdone!\n";
    print Log "########################################\n";
    print Log "blast Done and elapsed time : ",time()-$time,"s\n";
    print Log "########################################\n";
    $step++ ;
}

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
	
	eg:
	perl $Script -i -o -k -c

Usage:

    -ref	<file>	input ref files;must be give
    -out	<dir>	output dir;must be give
    -insert	<file>	insert fasta file
    -len    <number> input insert number;defult 400
    -fqlist <file>  fqlist file
    -step   pipiline control
    -stop   pipiline control
    -h			Help

USAGE
        print $usage;
        exit;
}
