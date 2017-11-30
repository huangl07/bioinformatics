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
		This script is write to change marker to bin file.
		Version: $ver

	Usage:

		-i           marker file                 <infile>     must be given
		-o           outfile prefix              <outfile>    must be given
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
my $markerfile = $opts{i} ;
my $outfile = $opts{o} ;

## reading marker file and get bin
my %hbin = ();
my %hbin_info = ();
my $asample_name = &parse_marker_file_and_get_bin($markerfile, \%hbin, \%hbin_info) ;

## out put
&out_put_bin_result($asample_name, \%hbin, \%hbin_info, $outfile);


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

#&parse_marker_file_and_get_bin($markerfile, \%hbin, \%hbin_info) ;
sub parse_marker_file_and_get_bin()
{
	my ($markerfile, $ahbin, $ahbin_info) = @_ ;

	open (IN, $markerfile) || die "Can't open $markerfile, $!\n" ;
	my @sample_name = ();
	my @pre_types = ();
	my $pre_chr = "" ;
	my $bin_num = 0 ;
	while (<IN>){
		chomp ;
		next if (m/^\s*$/) ;
		my ($marker_ID, $chr, $pos, $type, $P, $M, @sample_types) = split ;
		if ($marker_ID =~ /ID$/i){
			@sample_name = @sample_types ;
			next ;
		}
		if ($chr ne $pre_chr){ ## add new bin
			$bin_num++ ;
			my $start = $pos ;
			my $end = $pos ;
			my $marker_num = 1 ;
			push @{$ahbin_info->{$chr}}, [$bin_num, $start, $end, $marker_num, $marker_ID] ;
			push @{$ahbin->{$chr}}, \@sample_types ;
			$pre_chr = $chr ;
		}
		else{
			my $flag = &judge_bin_info($ahbin->{$chr}->[-1], \@sample_types);
			if ($flag == 1){
				$bin_num++ ;
				my $start = $pos ;
				my $end = $pos ;
				my $marker_num = 1 ;
				push @{$ahbin_info->{$chr}}, [$bin_num, $start, $end, $marker_num, $marker_ID] ;
				push @{$ahbin->{$chr}}, \@sample_types ;
			}
			else{
				$ahbin_info->{$chr}->[-1][2] = $pos ;
				$ahbin_info->{$chr}->[-1][3] += 1 ;
				push @{$ahbin_info->{$chr}->[-1]}, $marker_ID ;
			}
		}
	}
	close(IN);

	return(\@sample_name);
}

#my $flag = &judge_bin_info($ahbin->{$chr}->[-1], \@sample_types);
sub judge_bin_info()
{
	my ($abin, $asample_types) = @_ ;
	my $flag = 0 ;

	my @miss_index = ();
	for (my $i=0; $i<@{$abin}; $i++){
		if ($abin->[$i] ne $asample_types->[$i]){
			if ($abin->[$i] eq '--'){
				push @miss_index, $i ;
			}
			else{
				$flag = 1 ;
			}
		}
	}
	if ($flag == 0 && @miss_index >= 1){
		for my $index (@miss_index){
			$abin->[$index] = $asample_types->[$index] ;
		}
	}

	return($flag);
}

#&out_put_bin_result($asample_name, \%hbin, \%hbin_info, $outfile);
sub out_put_bin_result()
{
	my ($asample_name, $ahbin, $ahbin_info, $outfile) = @_ ;

	open (BIN, ">$outfile.bin") || die "Can't creat $outfile.bin, $!\n" ;
	print BIN "MarkerID\tChr\tpos\ttype\tP\tM\t", join("\t", @{$asample_name}), "\n" ;
	open (INF, ">$outfile.bininfo") || die "Can't creat $outfile.bininfo, $!\n" ;
	print INF "MarkerID\tChr\tpos\tmarker_num\tmarkers\n" ;
	for my $chr (sort keys %{$ahbin}){
		for (my $i=0; $i<@{$ahbin->{$chr}}; $i++){
			print BIN "Block",$ahbin_info->{$chr}->[$i][0], "\t$chr\t", $ahbin_info->{$chr}->[$i][1],"-", $ahbin_info->{$chr}->[$i][2] ;
			print BIN "\taaxbb\taa\tbb\t", join("\t", @{$ahbin->{$chr}->[$i]}), "\n" ;
			print INF "Block", $ahbin_info->{$chr}->[$i][0], "\t$chr\t", $ahbin_info->{$chr}->[$i][1],"-", $ahbin_info->{$chr}->[$i][2] ;
			print INF "\t", $ahbin_info->{$chr}->[$i][3], "\t", join(":", @{$ahbin_info->{$chr}->[$i]}[4..$#{$ahbin_info->{$chr}->[$i]}]),"\n" ;
		}
	}
	close(BIN);
	close(INF);

	return ;
}

#push @{$ahbin_info->{$chr}}, [$bin_num, $start, $end, $marker_num] ;

