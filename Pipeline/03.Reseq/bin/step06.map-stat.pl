#!/usr/bin/env perl 
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($proc,$vcflist,$bamlist,$chrlist,$annolist,$svlist,$cnvlist,$dOut,$dShell,$metric,$dictfile);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"bam:s"=>\$bamlist,
	"chr:s"=>\$chrlist,
	"proc:s"=>\$proc,
	"dict:s"=>\$dictfile,
	"out:s"=>\$dOut,
	"dsh:s"=>\$dShell,
	"met:s"=>\$metric,
			) or &USAGE;
&USAGE unless ($bamlist and $dOut and $dShell and $chrlist and $dictfile);
$proc||=20;
mkdir $dOut if (!-d $dOut);
mkdir $dShell if (!-d $dShell);
$bamlist=ABSOLUTE_DIR("$bamlist");
$dOut=ABSOLUTE_DIR("$dOut");
$dShell=ABSOLUTE_DIR("$dShell");
$chrlist=ABSOLUTE_DIR("$chrlist");
$metric=ABSOLUTE_DIR("$metric") if($metric ne "NULL");
$dictfile=ABSOLUTE_DIR("$dictfile");
open In,$metric;
my %metric;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sample,$metric)=split(/\t/,$_);
	$metric{$sample}=$metric;
}
close In;
open SH,">$dShell/06.map-stat.sh";
open In,$bamlist;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($sampleID,$bam)=split(/\t/,$_);
	print SH "samtools stats $bam > $dOut/$sampleID.all.mapstat && ";
	print SH "perl $Bin/bin/map_stat.pl -b $dOut/$sampleID.all.mapstat -m $metric{$sampleID} -i $dOut/$sampleID.insert -c $dOut/$sampleID.coverage -d $dictfile -o $dOut/$sampleID.result.stat -k $sampleID && ";
	print SH "Rscript $Bin/bin/insert_size.R --i $dOut/$sampleID.insert --o $dOut/$sampleID.insert && ";
	print SH "Rscript $Bin/bin/coverage_depth.R --i $dOut/$sampleID.coverage --o $dOut/$sampleID.depth \n";
	print SH "samtools depth $bam > $dOut/$sampleID.depth && ";
	print SH "$Bin/bin/depth_stat_windows -i $dOut/$sampleID.depth -o $dOut/$sampleID.depth.fordraw -w 100000 && ";
	print SH "Rscript $Bin/bin/genomeCoveragehorizontalArea.R --infile $dOut/$sampleID.depth.fordraw  --idfile $chrlist --outfile $dOut/$sampleID.genome.coverage --group.col 1 --x.col 2 --y.col 3 --x.lab Sequence-Position --y.lab AverageDepth-log2 --skip 0 --unit 100kb --log2\n";
}
close In;
close SH;
my $job="perl /mnt/ilustre/users/dna/.env/bin/qsub-sge.pl  --Resource mem=3G --CPU 1  --maxjob $proc $dShell/06.map-stat.sh";
print $job;
`$job`;
`paste $dOut/*.all.mapstat > $dOut/total.mapstat`;
my $who=`whoami`;
chomp($who);
SendAttachMail("mapping done!you must check $dOut/total.mapstat",'$dOut/total.mapstat',"who\@majorbio.com");

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
sub SendAttachMail
{
    my $mail_content = shift;
    my $mail_attach = shift;
	my $mail_to=shift;
	my $mail_from="long.huang\@majorbio.com";
	my $mail_cc="long.huang\@majorbio.com,dna\@majorbio.com";
	my $mail_subject="poor luck,qsub error";
    my $msg=MIME::Lite->new(
                    From=>$mail_from,
                    To=>$mail_to,
                    Cc=>$mail_cc,
                    Subject=>$mail_subject,
                    Type=>'TEXT',
                    Data=>$mail_content,);
	if ($mail_attach != "") {
		    $msg->attach(
            Type=>'AUTO',
            Path=>$mail_attach,
            Filename=>$mail_attach,);
	}
	$msg->send('smtp', "smtp.majorbio.com", AuthUser=>"long.huang\@majorbio.com", AuthPass=>"cwn.711hl" );
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
  -bam	<file>	input bam list 
  -chr	<file>	input chr list
  -proc <num>   number of process for qsub,default 20
  -met	<file>	input metric list
  -dict	<file>	input dict file
  -out	<dir>	output dir
  -dsh	<dir>	output shell dir
  -h         Help

USAGE
        print $usage;
        exit;
}
