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
my $Title="draw_haplotype";
my $version="v1.2.1";
my $queue="middle.q";
my $cpu=50;
my $vf="5G";
######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"test=s","i=s","od=s","p=s","g=s","n=s","sn=s","max=s","maxproc=s","f=s","t=s","d=s","b=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{od}) || !defined($opts{p}) || !defined($opts{g}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           marker or block genotype file            <infile>     must be given
		-od          result file dir                          <outdir>     must be given
		-p           outfile prefix                           <string>     must be given
		-g           grade of Population, if (g=2)            <int>        must be given
		             then type fo population if f2
		-b           bin info file, for extend bin file       [file]       optional
		-d           distance for one block for extend        [int]        optional [100][kb]
		-n           marker or block number draw in one dot   [int]        optional [1]
		-sn          split marker number in one file          [int]        optional [2000]
		-f           whether swap genotype, 1:yes,0:no        [int]        optional [1]
		-max         max chr/scaffold for draw                [int]        optional [500]
		-maxproc     max process for qsub                     [int]        optional [100]
		-t           merge marker or split marker             [int]        optional [1]
		             file, 1:split,0:merge
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
my $infile = $opts{i} ;
$infile = &ABSOLUTE_DIR($infile) ;
my $outdir = $opts{od} ;
mkdir $outdir if (!-d $outdir) ;
$outdir = &ABSOLUTE_DIR($outdir) ;
my $prefix = $opts{p} ;
my $grade = $opts{g} ;
my $bininfofile = defined $opts{b} ? $opts{b} : "" ;
if (defined $opts{b}){
	$bininfofile = &ABSOLUTE_DIR($bininfofile) ;
}
my $test=$opts{test};
createLog($Title,$version,$$,"$outdir/log/",$test);	



my $distance = defined $opts{d} ? $opts{d} : 100 ;
my $marker_number = defined $opts{n} ? $opts{n} : 1 ;
my $split_marker_number = defined $opts{sn} ? $opts{sn} : 2000 ;
my $max_chr_num = defined $opts{max} ? $opts{max} : 500 ;
my $maxproc = defined $opts{maxproc} ? $opts{maxproc} : 100 ;
my $swap_flag = defined $opts{f} ? $opts{f} : 1 ;
my $split_flag = defined $opts{t} ? $opts{t} : 1 ;

## reading marker or block file and split chr
stepStart(1,"reading marker or block file and split chr");
my $agenotype_files = &split_chr_genotype($infile, $outdir, $prefix, $max_chr_num, $swap_flag);
stepTime(1);
## extend bin file
stepStart(2,"extend bin file");
if (defined $opts{b}){
	$agenotype_files = &extend_bin_file($agenotype_files, $bininfofile, $distance, $outdir, $prefix);
}
stepTime(2);
## merge marker or split marker
stepStart(3,"merge marker or split marker");
my $amerge_genotype_files ;
if ($split_flag == 0){
	$amerge_genotype_files = &merge_marker($agenotype_files, $outdir, $prefix, $marker_number);
}
else{
	$amerge_genotype_files = &split_marker($agenotype_files, $outdir, $prefix, $split_marker_number);
}
stepTime(3);
## get loc files
stepStart(4,"get loc files");
my $aloc_files = &get_loc_file($amerge_genotype_files, $outdir, $prefix, $grade);
stepTime(4);
## draw haplotype map
stepStart(5,"draw haplotype map");
&draw_haplotype_map($aloc_files, $outdir, $prefix);
stepTime(5);

totalTime();	
###############Time
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

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

#my $agenotype_files = &split_chr_genotype($infile, $outdir, $prefix, $max_chr_num, $swap_flag);
sub split_chr_genotype()
{
	my ($infile, $outdir, $prefix, $max_chr_num, $swap_flag) = @_ ;

	my $dir = "$outdir/chr_genotype" ;
	mkdir $dir if (!-d $dir) ;
	my $cmd = "perl $Bin/get_chr_genotype.pl -g $infile -od $dir -p $prefix -max $max_chr_num -f $swap_flag" ;
	runOrDie($cmd);
	#&run_or_die($cmd) ;

	my @genotype_files = glob("$dir/$prefix.*.genotype");

	return(\@genotype_files);
}

#$agenotype_files = &extend_bin_file($agenotype_files, $bininfofile, $distance, $outdir, $prefix);
sub extend_bin_file()
{
	my ($agenotype_files, $bininfofile, $distance, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/extend_binfile" ;
	mkdir $dir if (!-d $dir) ;
	my $job_num = 0 ;
	my @extend_bin_files = ();

	my $shfile = "$dir/extend_binfile.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $genotypefile (@{$agenotype_files}){
		my $basename = basename($genotypefile);
		my $outfile = "$dir/$basename" ;
		my $cmd = "perl $Bin/extend_bin.pl -i $genotypefile -f $bininfofile -o $outfile -d $distance" ;
		print SH $cmd, "\n" ;
		$job_num++ ;
		push @extend_bin_files, "$outfile.extend" ;
	}
	close(SH) ;
	qsubOrDie($shfile,$queue,$cpu,$vf);	
	#&qsub_mult($shfile, $maxproc, $job_num);

	return(\@extend_bin_files);
}

#my $merge_genotype_files = &merge_marker($agenotype_files, $outdir, $prefix, $marker_number);
sub merge_marker()
{
	my ($agenotype_files, $outdir, $prefix, $marker_number) = @_ ;

	my $dir = "$outdir/merge_marker" ;
	mkdir $dir if (!-d $dir) ;
	my $job_num = 0 ;
	my @merge_genotype_files = ();

	my $shfile = "$dir/merge_marker.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $genotypefile (@{$agenotype_files}){
		my $basename = basename($genotypefile) ;
		$basename =~ s/.genotype$/.merge.genotype/ ;
		my $outfile = "$dir/$basename" ;
		my $cmd = "perl $Bin/merge_marker.pl -i $genotypefile -o $outfile -n $marker_number" ;
		print SH $cmd,"\n" ;
		$job_num++ ;
		push @merge_genotype_files, $outfile ;
	}
	close(SH);
	qsubOrDie($shfile,$queue,$cpu,$vf);
	#&qsub_mult($shfile, $maxproc, $job_num);

	return (\@merge_genotype_files);
}

#$amerge_genotype_files = &split_marker($agenotype_files, $outdir, $prefix, $split_marker_number);
sub split_marker()
{
	my ($agenotype_files, $outdir, $prefix, $split_marker_number) = @_ ;

	my $dir = "$outdir/split_marker" ;
	mkdir $dir if (!-d $dir) ;
	my $job_num = 0 ;

	my $shfile = "$dir/split_marker.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $genotypefile (@{$agenotype_files}){
		my $basename = basename($genotypefile) ;
		my $new_dir = "$dir/$basename" ;
		mkdir $new_dir if (!-d $new_dir) ;
		my $cmd = "perl $Bin/split_marker.pl -i $genotypefile -od $new_dir -p $basename -n $split_marker_number" ;
		print SH $cmd, "\n" ;
		$job_num++ ;
	}
	close(SH);
	qsubOrDie($shfile,$queue,$cpu,$vf);
#	&qsub_mult($shfile, $maxproc, $job_num);

	my @split_genotype_file = glob("$dir/*/*.split*");

	return(\@split_genotype_file);
}

#my $loc_files = &get_loc_file($amerge_genotype_files, $outdir, $prefix, $grade);
sub get_loc_file()
{
	my ($amerge_genotype_files, $outdir, $prefix, $grade) = @_ ;

	my $dir = "$outdir/loc_file" ;
	mkdir $dir if (!-d $dir) ;
	my $job_num = 0 ;
	my @loc_files = ();

	my $shfile = "$dir/loc_file.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $infile (@{$amerge_genotype_files}){
		my $basename = basename($infile);
		my $outfile = "$dir/$basename" ;
		my $cmd = "perl $Bin/genotype2loc_Rils_v1.2.pl -i $infile -o $outfile -g $grade -add" ;
		print SH $cmd, "\n" ;
		push @loc_files, "$outfile.loc" ;
		$job_num++ ;
	}
	close(SH);
	qsubOrDie($shfile,$queue,$cpu,$vf);
	#&qsub_mult($shfile, $maxproc, $job_num);
	
	return(\@loc_files);
}

#&draw_haplotype_map($aloc_files, $outdir, $prefix);
sub draw_haplotype_map()
{
	my ($aloc_files, $outdir, $prefix) = @_ ;

	my $dir = "$outdir/haplotype_map" ;
	mkdir $dir if (!-d $dir);
	my $job_num = 0 ;

	my $shfile = "$dir/haplotype_map.sh" ;
	open (SH, ">$shfile") || die "Can't creat $shfile, $!\n" ;
	for my $loc_file (@{$aloc_files}){
		my $basename = basename($loc_file);
		my $outfile = "$dir/$basename" ;
		my $cmd = "perl $Bin/draw_haplotype_map.pl -i $loc_file -o $outfile" ;
		print SH $cmd, "\n" ;
		$job_num++ ;
	}
	close(SH);
	qsubOrDie($shfile,$queue,$cpu,$vf);
	#&qsub_mult($shfile, $maxproc, $job_num);

	return ;
}


