#!/usr/bin/perl -w
# 
# Copyright (c) BIO_MARK 2014
# Writer:         Dengdj <dengdj@biomarker.com.cn>
# Program Date:   2014
# Modifier:       Dengdj <dengdj@biomarker.com.cn>
# Last Modified:  2014.
my $ver="1.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use newPerlBase;
my $Title="genetic_map:rseq_genticmap";	
my $version="v1.2.1";	
######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"test=s","i=s","t=s","od=s","k=s","cn=s","g=s","P=s","M=s","Pd=s","Md=s","mw=s","w=s","s=s","ad=s","bd=s","n=s","fl=s","imp","cor","bin","h","linkage" );
	
#&help()if(defined $opts{h});
if((!defined($opts{i}) && !defined($opts{od})) || !defined($opts{od}) || !defined($opts{k}) || !defined($opts{cn}) || !defined($opts{g}) || !defined($opts{P}) || !defined($opts{M}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver
		v1.1	use new imputation or correction method.

	Usage:

		-i           vcf file                       <infile>     must have one with -t
		-t           genotype file                  <infile>     must have one with -i
		-od          result file dir                <outdir>     must be given
		-k           outfile prefix                 <string>     must be given
		-cn          chr number                     <int>        must be given
		-g           grade of Population            <int>        must be given
		-P           paternal sample ID             <string>     must be given
		-M           maternal sample ID             <string>     must be given
		-Pd          paternal min depth             [int]        optional [4]
		-Md          maternal min depth             [int]        optional [4]
		-mw          match word for filter          [string]     optional ['Chr']
		-w           window length for imputation   [int]        optional [15]
		-s           step length for imputation     [int]        optional [1]
		-ad          aa depth                       [int]        optional [11]
		-bd          bb depth                       [int]        optional [11]
		-n           marker number for correct      [int]        optional [50]
		-fl          filter len for bin             [int]        optional [100][kb]
		-imp         do imputation                  [null]       optional
		-cor         do correction                  [null]       optional
		-bin         do bin fined and filter        [null]       optional
		-linkage     do you want linkage                         optional
		-h           Help document

	Usage End.

	exit;
}

###############Time
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
################
## get parameters
my $vcffile = defined $opts{i} ? $opts{i} : "" ;
if (defined $opts{i}){
	$vcffile = &ABSOLUTE_DIR($vcffile) ;
}
my $genotypefile = defined $opts{t} ? $opts{t} : "" ;
if (defined $opts{t}){
	$genotypefile = &ABSOLUTE_DIR($genotypefile);
}
my $outdir = $opts{od} ;
mkdir $outdir if (!-d $outdir) ;
$outdir = &ABSOLUTE_DIR($outdir) ;
my $prefix = $opts{k} ;
my $chr_num = $opts{cn} ;
my $grade = $opts{g} ;
my $paternal_ID = $opts{P} ;
my $maternal_ID = $opts{M} ;
my $paternal_depth = defined $opts{Pd} ? $opts{Pd} : 4 ;
my $maternal_depth = defined $opts{Md} ? $opts{Md} : 4 ;
my $match_word = defined $opts{mw} ? $opts{mw} : 'Chr' ;
my $window_len = defined $opts{w} ? $opts{w} : 15 ;
my $step_len = defined $opts{s} ? $opts{s} : 1 ;
my $a_depth = defined $opts{ad} ? $opts{ad} : 11 ;
my $b_depth = defined $opts{bd} ? $opts{bd} : 11 ;
my $n = defined $opts{n} ? $opts{n} : 50 ;
my $filter_len = defined $opts{fl} ? $opts{fl} : 100 ;
my $test=$opts{test};
createLog($Title,$version,$$,"$outdir/log/",$test);

## vcf2marker
stepStart(1,"vcf2marker"); 
my $raw_marker_file = "" ;
if (defined $opts{i}){
	$raw_marker_file = &vcf2marker($vcffile, $outdir, $prefix);
}
else{
	$raw_marker_file = $genotypefile ;
}
stepTime(1);
## imputation
stepStart(2,"imputation"); 
my $imp_marker_file = $raw_marker_file ;
if (defined $opts{imp}){
	$imp_marker_file = &imputation_marker($raw_marker_file, $outdir, $prefix);
}
stepTime(2);
## correct
stepStart(3,"correct"); 
my $corr_marker_file = $imp_marker_file ;
if (defined $opts{cor}){
	$corr_marker_file = &corr_block($imp_marker_file, $outdir, $prefix);
}
stepTime(3);
## get bin and filter
stepStart(4,"get bin and filter");
my $binfile = $corr_marker_file ;
if (defined $opts{bin}){
	$binfile = &get_marker_bin($corr_marker_file, $outdir, $prefix);
}
stepTime(4);
## construct genetic map
if(defined $opts{"linkage"}){
	stepStart(5," construct genetic map");
	&construct_genetic_map($binfile, $outdir, $prefix);
	stepTime(5);
}
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";
totalTime();	
###############Subs
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;
	$cur_dir =~ s/\n$//;
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;
		$dir =~ s/\n$// ;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;
		$return =~ s/\n$// ;
	}
	else
	{
		warn "Warning just for file and dir [$in]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# &show_log("txt")
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}

#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	my $start_time = &show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		my $end_time = &show_log("Error: command fail: $cmd");
		&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
		exit(1);
	}
	my $end_time = &show_log("done.");
	&show_log("Elaseped time: ".($end_time - $start_time)."s\n");
	return ;
}

## qsub
sub qsub()
{
	my ($shfile, $maxproc) = @_ ;
	if (`hostname` =~ /cluster/){
		my $cmd = "qsub-sge.pl --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}
	else{
		my $cmd = "ssh cluster -Y qsub-sge.pl --maxproc $maxproc --reqsub $shfile --independent" ;
		&run_or_die($cmd);
	}

	return ;
}

## qsub_mult($shfile, $max_proc, $job_num)
sub qsub_mult()
{
	my ($shfile, $max_proc, $job_num) = @_ ;
	if ($job_num > 500){
		my @shfiles = &cut_shfile($shfile);
		for my $file (@shfiles){
			&qsub($file, $max_proc);
		}
	}
	else{
		&qsub($shfile, $max_proc) ;
	}
}

#my @shfiles = &cut_shfile($shfile);
sub cut_shfile()
{
	my ($file) = @_ ;
	my @files = ();
	my $num = 0 ;
	open (IN, $file) || die "Can't open $file, $!\n" ;
	(my $outprefix = $file) =~ s/.sh$// ;
	while (<IN>){
		chomp ;
		if ($num % 500 == 0){
			close(OUT);
			my $outfile = "$outprefix.sub_".(int($num/500)+1).".sh" ;
			push @files, $outfile ;
			open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
		}
		print OUT $_, "\n" ;
		$num ++ ;
	}
	close(IN);
	close(OUT);

	return @files ;
}

#my $raw_file = &vcf2marker($vcffile, $outdir, $prefix);
sub vcf2marker()
{
	my ($vcffile, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/marker_info" ;
	mkdir $dir if (!-d $dir) ;

	# vcf 2 marker
	my $marker_file = "$dir/$prefix.genotype" ;
	my $cmd = "perl $Bin/vcf2marker.pl -i $vcffile -o $marker_file -P $paternal_ID -M $maternal_ID -Pd $paternal_depth -Md $maternal_depth" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	# filter scaffold
	my $outfile = "$dir/$prefix.Chr.genotype" ;
	$cmd = "perl $Bin/retain_chr.pl -i $marker_file -o $outfile -m $match_word" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	# stat genotype
	&stat_genotype($outfile);

	return($outfile);
}

#my ($imp_marker_file) = &imputation_marker($raw_marker_file, $outdir, $prefix);
sub imputation_marker()
{
	my ($raw_marker_file, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/imputation" ;
	mkdir $dir if (!-d $dir) ;

	#my $basename = basename($raw_marker_file) ;
	#$basename =~ s/.genotype$/.imp.genotype/ ;
	my $imp_marker_file = "$dir/$prefix.imp" ;
	my $cmd = "perl $Bin/marker_imputation_v1.1.pl -i $raw_marker_file -od $dir -p $prefix -w $window_len -s $step_len -ad $a_depth -bd $b_depth" ;
	runOrDie($cmd);
	#my $cmd = "perl $Bin/map_corr.pl -i $raw_marker_file -o $imp_marker_file -w $window_len -s $step_len" ;
	#&run_or_die($cmd);

	# stat genotype
	&stat_genotype($imp_marker_file);

	# draw haplotype map
	#&draw_haplotype_map($imp_marker_file, "$dir/haplotype_map", $prefix);

	return($imp_marker_file);
}

#my $corr_block_file = &corr_block($raw_block_file, $outdir, $prefix)
sub corr_block()
{
	my ($imp_marker_file, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/correction" ;
	mkdir $dir if (!-d $dir) ;

	my $basename = basename($imp_marker_file) ;
	my $corr_marker_file = "$dir/$basename.corr" ;
	my $cmd = "perl $Bin/marker_correct.pl -i $imp_marker_file -o $corr_marker_file -n $n" ;
	#my $cmd = "perl $Bin/block_corr.pl -i $raw_block_file -o $corr_block_file" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	return($corr_marker_file);
}

#$binfile = &get_marker_bin($corr_marker_file, $outdir, $prefix);
sub get_marker_bin()
{
	my ($marker_file, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/marker_binfile" ;
	mkdir $dir if (!-d $dir);

	my $basename = basename($marker_file);
	my $binfile = "$dir/$basename" ;
	my $cmd = "perl $Bin/marker2bin.pl -i $marker_file -o $binfile" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	my $filter_binfile = "$dir/$basename.filter" ;
	$cmd = "perl $Bin/filter_bin.pl -i $binfile.bin -f $binfile.bininfo -o $filter_binfile -l $filter_len" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	return("$filter_binfile.bin");
}

#&construct_genetic_map($corr_block_file, $outdir, $prefix);
sub construct_genetic_map()
{
	my ($corr_block_file, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/genetic_map" ;
	mkdir $dir if (!-d $dir) ;

	# get position file
	my $posi_file = "$dir/$prefix.posi" ;
	my $cmd = "perl $Bin/marker2posi.pl -i $corr_block_file -o $posi_file" ;
	runOrDie($cmd);
	#&run_or_die($cmd);

	# change block file
	my $block_file = "$dir/$prefix.block" ;
	$cmd = "perl $Bin/change_block_format.pl -i $corr_block_file -o $block_file" ;
	runOrDie($cmd);
#	&run_or_die($cmd);

	# genetic map
	my $logfile = "$dir/$prefix.log" ;
##------------------------------------------------update the linkage_inbred by liangsh
	$cmd = "perl $Bin/Linkage/linkage_Inbred/linkage_Inbred.pl -i $block_file -k $prefix -d $dir -ref $posi_file -popt Ri$grade -nChro $chr_num -MSTmap > $logfile" ;
	runOrDie($cmd);
	#&run_or_die($cmd);

	return ;
}

#&stat_genotype($outfile);
sub stat_genotype()
{
	my ($infile) = @_ ;
	my $outfile = "$infile.ratio.stat" ;
	my $cmd = "perl $Bin/stat/marker_ratio_stat.pl -i $infile -o $outfile" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	return ;
}

#&draw_haplotype_map($imp_marker_file, "$dir/haplotype_map", $prefix);
sub draw_haplotype_map()
{
	my ($marker_file, $dir) = @_ ;

	my $cmd = "perl $Bin/draw_haplotype/draw_haplotype_map_process.pl -i $marker_file -od $dir -p $prefix -g $grade" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	return ;
}

