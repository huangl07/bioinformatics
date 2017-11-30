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

######################请在写程序之前，一定写明时间、程序用途、参数说明；每次修改程序时，也请做好注释工作

my %opts;
GetOptions(\%opts,"i=s","f=s","o=s","l=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{f}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to filter bin by bin length or segregation distortion.
		Version: $ver
		v1.0:   Filter bin by length.

	Usage:
		-i           bin file                 <infile>     must be given
		-f           bin info file            <infile>     must be given
		-o           outfile prefix           <outfile>    must be given
		-l           length for filter        [int]        optional [100][kb]
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
my $binfile = $opts{i} ;
my $infofile = $opts{f} ;
my $outfile = $opts{o} ;
my $min_len = defined $opts{l} ? $opts{l}*1000 : 100*1000 ;

## reading bininfo file and do filter
my %hinfo = ();
&filter_bininfofile_by_length($infofile, \%hinfo, $outfile);

## filter bin file
&filter_binfile($binfile, \%hinfo, $outfile);

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

#&filter_bininfofile_by_length($infofile, \%hinfo, $outfile);
sub filter_bininfofile_by_length()
{
	my ($infofile, $ahinfo, $outfile) = @_ ;

	open (IN, $infofile) || die "Can't open $infofile, $!\n" ;
	open (OUT, ">$outfile.bininfo") || die "Can't creat $outfile.bininfo, $!\n" ;
	my $filter_num = 0 ;
	my $total_num = 0 ;
	while (<IN>){
		chomp ;
		next if(m/^\s*$/ || m/^\#/);
		my ($markerID, $chr, $pos, $marker_num, $markers) = split ;
		if ($markerID =~ m/MarkerID/i || $markerID =~ m/BlockID/i){
			print OUT $_,"\n" ;
			next ;
		}
		$total_num++ ;
		my ($start, $end) = split /-/, $pos ;
		if ($end - $start < $min_len){
			$filter_num++ ;
			next ;
		}
		$ahinfo->{$markerID} = 1 ;
		print OUT $_,"\n" ;
	}
	close(IN);
	close(OUT);
	open (FT, ">$outfile.stat") || die "Can't creat $outfile.log, $!\n" ;
	print FT "Total:\t$total_num\n" ;
	print FT "Filter:\t$filter_num\n" ;
	close(FT);

	return ;
}

#&filter_binfile($binfile, \%hinfo, $outfile);
sub filter_binfile()
{
	my ($binfile, $ahinfo, $outfile) = @_ ;

	open (IN, $binfile) || die "Can't open $binfile, $!\n" ;
	open (OUT, ">$outfile.bin") || die "Can't creat $outfile.bin, $!\n" ;
	while (<IN>){
		chomp ;
		next if(m/^\s*$/ || m/^\#/);
		my ($markerID, $chr, $pos, $marker_num, $markers) = split ;
		if ($markerID =~ m/MarkerID/i || $markerID =~ m/BlockID/i){
			print OUT $_,"\n" ;
			next ;
		}
		if (defined $ahinfo->{$markerID}){
			print OUT $_,"\n" ;
		}
	}
	close(IN);
	close(OUT);

	return ;
}


