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
GetOptions(\%opts,"i=s","f=s","o=s","d=s","h");

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{f}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to extend bin file.
		Version: $ver

	Usage:
		-i           bin file                         <infile>     must be given
		-f           bin info file                    <infile>     must be given
		-o           out prefix                       <outfile>    must be given
		-d           distance for one block           [int]        optional [100][kb]
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
my $distance = defined $opts{d} ? $opts{d}*1000 : 100*1000 ;

## reading bin info file and get bin extend number
my %hextend = ();
&get_bin_extend_num($infofile, \%hextend, $outfile);

## reading bin file and extend
&extend_bin(\%hextend, $binfile, $outfile);

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

#&get_bin_extend_num($infofile, \%hextend, $outfile);
sub get_bin_extend_num()
{
	my ($infofile, $ahextend, $outfile) = @_ ;

	open (IN, $infofile) || die "Can't open $infofile, $!\n" ;
	open (MID, ">$outfile.mid") || die "Can't creat $outfile.mid, $!\n" ;
	print MID "MarkerID\tChr\tpos\tmarker_num\textend_num\n" ;
	while (<IN>){
		next if (m/^\#/ || m/^\s*$/);
		next if (m/MarkerID/i) ;
		my ($markerID, $chr, $pos, $marker_num) = split ;
		my ($start, $end) = split /-/, $pos ;
		my $num = int(($end - $start)/$distance)+1 ;
		$ahextend->{$markerID} = $num ;
		print MID "$markerID\t$chr\t$pos\t$marker_num\t$num\n" ;
	}
	close(IN);
	close(MID);

	return ;
}

#&extend_bin(\%hextend, $binfile, $outfile);
sub extend_bin()
{
	my ($ahextend, $binfile, $outfile) = @_ ;

	open (IN, $binfile) || die "Can't open $binfile, $!\n" ;
	open (OUT, ">$outfile.extend") || die "Can't creat $outfile.extend, $!\n" ;
	while(<IN>){
		chomp ;
		next if (m/^\*$/);
		(my $markerID) = (split)[0] ;
		if (defined $ahextend->{$markerID}){
			print OUT "$_\n"x$ahextend->{$markerID} ;
		}
		else{
			print OUT "$_\n" ;
		}
	}
	close(IN) ;
	close(OUT) ;

	return ;
}


