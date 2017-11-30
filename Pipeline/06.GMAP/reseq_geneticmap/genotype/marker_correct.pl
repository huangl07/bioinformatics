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
GetOptions(\%opts,"i=s","o=s","n=s","h");

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to correct imputation marker file.
		Version: $ver
		v1.0: Function realization. Do not consider '--' specially.

	Usage:
		-i           marker file                           <infile>     must be given
		-o           corrected file                        <outfile>    must be given
		-n           min marker number                     <int>        must be given
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
my $marker_file = $opts{i} ;
my $outfile = $opts{o} ;
my $min_num = $opts{n} ;

## reading marker file
my %hmarker = ();
my %hmarker_info = ();
my %hmarker_num = ();
my $asample_name = &reading_marker_file($marker_file, \%hmarker, \%hmarker_info, \%hmarker_num);

## correct marker
my %hwindow = ();
&corrcet_bin_marker(\%hmarker, \%hmarker_info, \%hwindow);

## out put marker
&output($asample_name, \%hmarker_info, \%hwindow, $outfile) ;


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

#&reading_marker_file($marker_file, \%hmarker, \%hmarker_info, \%hmarker_num);
sub reading_marker_file()
{
	my ($marker_file, $ahmarker, $ahmarker_info, $ahmarker_num) = @_ ;

	open (IN, $marker_file) || die "Can't open $marker_file, $!\n" ;
	my @sample_name = ();
	while (<IN>){
		chomp ;
		next if (m/^\s*$/);
		my ($marker_ID, $chr, $pos, $type, $P, $M, @sample_types) = split ;
		if ($marker_ID =~ /ID$/i){
			@sample_name = @sample_types ;
			next ;
		}
		for (my $i=0; $i<@sample_types; $i++){
			my $type = &swap_type($sample_types[$i]);
			push @{$ahmarker->{$sample_name[$i]}{$chr}}, $type ;
		}
		push @{$ahmarker_info->{$chr}}, [$marker_ID, $pos];
		$ahmarker_num->{$chr} ++ ;
	}
	close(IN);

	return(\@sample_name);
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

#&corrcet_bin_marker(\%hmarker, \%hmarker_info, \%hwindow);
sub corrcet_bin_marker()
{
	my ($ahmarker, $ahmarker_info, $ahwindow) = @_ ;

	for my $sample (keys %{$ahmarker}){
		#my $num = 0 ;
		#my $corr_num = 0 ;
		for my $chr (keys %{$ahmarker->{$sample}}){
			my $awin_marker = &get_win_marker($ahmarker->{$sample}->{$chr}) ;
			#$num += &marker_count($awin_marker);
			my $acorr_marker = &corr_win_marker($awin_marker) ;
			#$corr_num += &marker_count($acorr_marker);
			$ahwindow->{$sample}{$chr} = $acorr_marker ;
		}
		#print "$sample\t$num\t$corr_num\n" ;
	}

	return ;
}

#$num = &marker_count($awin_marker);
sub marker_count()
{
	my ($amarker) = @_ ;
	my $num = 0 ;

	for my $marker_info (@{$amarker}){
		$num += $marker_info->[1] ;
	}

	return ($num);
}


#&get_win_marker($ahmarker->{$sample}->{$chr}, $chr, $ahwindow) ;
sub get_win_marker()
{
	my ($amarker) = @_ ;

	my @win_marker = ();
	my $pre_marker = $amarker->[0] ;
	my $pre_index = 0 ;
	for (my $i=1; $i<@{$amarker}; $i++){
		my $marker = $amarker->[$i] ;
		if ($marker ne $pre_marker){
			push @win_marker, [$pre_marker, $i-$pre_index];
			$pre_marker = $marker ;
			$pre_index = $i ;
		}
	}
	push @win_marker, [$pre_marker, @{$amarker}-$pre_index] ;

	return (\@win_marker);
}

#my $corr_marker = &corr_win_marker($win_marker) ;
sub corr_win_marker()
{
	my ($awin_marker) = @_ ;
	my @corr_marker = ();

	if (@{$awin_marker} == 1){
		push @corr_marker, [$awin_marker->[0][0], $awin_marker->[0][1]] ;
		return (\@corr_marker) ;
	}

	my $start_index = 1 ;
	if ($awin_marker->[0][1] < $min_num && $awin_marker->[1][1] > $min_num){
		push @corr_marker, [$awin_marker->[1][0], $awin_marker->[0][1]+$awin_marker->[1][1]];
		$start_index += 1 ;
	}
	else{
		push @corr_marker, [$awin_marker->[0][0], $awin_marker->[0][1]] ;
	}
	for (my $i=$start_index; $i<@{$awin_marker}-1; $i++){
		my ($marker, $num) = @{$awin_marker->[$i]} ;
		if ($marker eq $corr_marker[-1][0]){
			$corr_marker[-1][1] += $awin_marker->[$i][1] ;
			next ;
		}
		if ($num < $min_num){
			if ($corr_marker[-1][0] eq $awin_marker->[$i+1][0]){
				$corr_marker[-1][1] += $awin_marker->[$i][1] ;
			}
			else{
				if ($marker eq 'ab'){
					push @corr_marker, [$awin_marker->[$i][0], $awin_marker->[$i][1]];
				}
				else{
					if ($corr_marker[-1][1] > $awin_marker->[$i+1][1]){
						$corr_marker[-1][1] += $num ;
					}
					else{
						$awin_marker->[$i+1][1] += $num ;
					}
				}
			}
		}
		else{
			push @corr_marker, [$awin_marker->[$i][0], $awin_marker->[$i][1]];
		}
	}
	if ($awin_marker->[-1][1] < $min_num){
		$corr_marker[-1][1] += $awin_marker->[-1][1] ;
	}
	elsif($corr_marker[-1][0] eq $awin_marker->[-1][0]){
		$corr_marker[-1][1] += $awin_marker->[-1][1] ;
	}
	else{
		push @corr_marker, [$awin_marker->[-1][0], $awin_marker->[-1][1]];
	}

	return(\@corr_marker);
}

#&output($asample_name, \%hmarker_info, \%hwindow, $outfile) ;
sub output()
{
	my ($asample_name, $ahmarker_info, $ahwindow, $outfile) = @_ ;

	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "MarkerID\tChr\tpos\ttype\tP\tM\t", join("\t", @{$asample_name}), "\n" ;
	#for my $sample (@{$asample_name}){
	#	my $num = 0 ;
	#	for my $chr (keys %{$ahwindow->{$sample}}){
	#		$num += &marker_count($ahwindow->{$sample}{$chr});
	#	}
	#	print "$sample\t$num\n" ;
	#}
	for my $chr (sort keys %{$ahmarker_info}){
		my @marker_info = ();
		for my $sample (@{$asample_name}){
			my $index = 0 ;
			for my $awin (@{$ahwindow->{$sample}{$chr}}){
				my ($marker, $num) = @{$awin} ;
				for (my $i=0; $i<$num; $i++){
					$marker_info[$index] .= "\t$marker" ;
					$index++ ;
				}
			}
			#print "$sample\t$index\n" ;
		}
		for (my $i=0; $i<@marker_info; $i++){
			print OUT $ahmarker_info->{$chr}[$i][0], "\t$chr\t", $ahmarker_info->{$chr}[$i][1],"\taaxbb\taa\tbb", $marker_info[$i],"\n" ;
		}
	}
	close(OUT);

	return ;
}


