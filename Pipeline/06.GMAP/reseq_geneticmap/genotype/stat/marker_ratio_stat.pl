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
GetOptions(\%opts,"i=s","o=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to statistic ratio of aa:ab:bb marker
		Version: $ver

	Usage:
		-i           marker file           <infile>     must be given
		-o           stat result file      <outfile>    must be given
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
my $infile  = $opts{i} ;
my $outfile = $opts{o} ;

## reading marker file and stat
open (IN,"$infile") || die "Can't open $infile,$!\n" ;
my @sample_name = ();
my %htype = ();
my $total_marker = 0 ;
while (<IN>){
	chomp ;
	next if (m/^\s*$/) ;
	my ($marker_id, $chr, $pos, $type, $P, $M, @sample_types) = split ;
	if ($marker_id =~ m/ID$/i){
		@sample_name = @sample_types ;
		next ;
	}
	for (my $i=0; $i<@sample_types; $i++){
		my $type = &swap_type($sample_types[$i]);
		$htype{$sample_name[$i]}{$type} ++ ;
	}
	$total_marker++ ;
}
close(IN);
open (OUT,">$outfile") || die "Can't creat $outfile,$!\n" ;
print OUT "#Sample\tTotal\tMis_num\tMis(%)\taa_num\tab_num\tbb_num\taa(%)\tab(%)\tbb(%)\n" ;
for (my $i=0; $i<@sample_name; $i++){
	my $name = $sample_name[$i] ;
	my $aa_marker = defined $htype{$name}{aa} ? $htype{$name}{aa} : 0 ;
	my $ab_marker = defined $htype{$name}{ab} ? $htype{$name}{ab} : 0 ;
	my $bb_marker = defined $htype{$name}{bb} ? $htype{$name}{bb} : 0 ;
	my $mis_marker = defined $htype{$name}{'--'} ? $htype{$name}{'--'} : 0 ;
	my $mis_percent = (int(($mis_marker/$total_marker)*10000))/100 ;
	my $persent_marker = $aa_marker + $ab_marker + $bb_marker ;
	my $aa_percent = (int(($aa_marker/$persent_marker)*10000))/100 ;
	my $ab_percent = (int(($ab_marker/$persent_marker)*10000))/100 ;
	my $bb_percent = (int(($bb_marker/$persent_marker)*10000))/100 ;
	print OUT "$name\t$total_marker\t$mis_marker\t$mis_percent\t$aa_marker\t$ab_marker\t$bb_marker\t$aa_percent\t$ab_percent\t$bb_percent\n" ;
}
close(OUT);
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

#my $type = &swap_type($sample_types[$i]);
sub swap_type()
{
	my ($raw_type) = @_ ;

	if ($raw_type =~ /-/){
		return ('--');
	}
	elsif ($raw_type =~ /\,/){
		my ($t1, $d1, $t2, $d2) = split /\,/, $raw_type ;
		if (!defined $t2){
			return ($t1.$t1);
		}
		else{
			return ($t1.$t2);
		}
	}
	else{
		return ($raw_type);
	}
}

