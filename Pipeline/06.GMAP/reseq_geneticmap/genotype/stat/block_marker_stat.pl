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
GetOptions(\%opts,"i=s","o=s","m=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{m}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver

	Usage:

		-i           posi file of block            <infile>     must be given
		-m           marker file of SNP            <infile>     must be given
		-o           result                        <outfile>    must be given
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
my $posifile  = $opts{i} ;
my $markerfile = $opts{m} ;
my $outfile = $opts{o} ;

## reading posi file
my %hposi = ();
my %hbin = ();
&reading_posi_file($posifile, \%hposi, \%hbin);

## reading marker file and stat
my %hstat = ();
&reading_marker_file($markerfile, \%hposi, \%hstat);

## out put
&output_result($outfile, \%hstat, \%hbin);

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

#&reading_posi_file($posifile, \%hposi, \%hbin);
sub reading_posi_file()
{
	my ($posifile, $ahposi, $ahbin) = @_ ;

	open (IN, $posifile) || die "Can't open $posifile, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\s$/ || m/^\#/);
		my ($id, $chr, $start, $end) = split ;
		next if ($id =~ /ID/i);
		push @{$ahposi->{$chr}}, [$start, $end] ;
		$ahbin->{$chr}++ ;
	}
	close(IN);

	return ;
}

#&reading_marker_file($markerfile, \%hposi, \%hstat);
sub reading_marker_file()
{
	my ($markerfile, $ahposi, $ahstat) = @_ ;

	my %hindex = ();
	for my $chr (keys %{$ahposi}){
		$hindex{$chr} = 0 ;
	}

	open (IN, $markerfile) || die "Can't open $markerfile, $!\n" ;
	while (<IN>){
		chomp ;
		next if (m/^\#/ || m/^\s*$/);
		my ($id, $chr, $pos) = split ;
		next if ($id =~ m/ID/i);
		for (my $i=$hindex{$chr}; $i<@{$ahposi->{$chr}}; $i++){
			my ($start, $end) = @{$ahposi->{$chr}->[$i]} ;
			if ($pos < $start){
				last ;
			}
			elsif ($pos <= $end){
				$ahstat->{$chr}++ ;
				$hindex{$chr} = $i ;
				last ;
			}
			else{
				$hindex{$chr} = $i ;
				next ;
			}
		}
	}
	close(IN);

	return ;
}

#&output_result($outfile, \%hstat, \%hbin);
sub output_result()
{
	my ($outfile, $ahstat, $ahbin) = @_ ;

	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "#Chr\tMarkerNum\tBinNum\n" ;
	for my $chr (sort keys %{$ahstat}){
		print OUT $chr,"\t", $ahstat->{$chr},"\t", $ahbin->{$chr},"\n" ;
	}
	close(OUT) ;

	return ;
}


