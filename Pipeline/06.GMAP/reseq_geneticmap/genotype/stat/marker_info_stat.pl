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
		This script is write to statistic marker info of each chr for each sample.
		Version: $ver

	Usage:
		-i           marker file            <infile>     must be given
		-o           stat file              <outfile>    must be given
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
my $outfile = $opts{o} ;

## reading infile and store each type for sample
my %htype = ();
my $asample_name = &get_marker_info($infile, \%htype);

## stat
&stat_marker_info($outfile, \%htype, $asample_name);

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

#&get_marker_info($infile, \%htype);
sub get_marker_info()
{
	my ($infile, $ahtype) = @_ ;

	open (IN, $infile) || die "Can't open $infile, $!\n" ;
	my @sample_name = ();
	while (<IN>){
		chomp ;
		next if (m/^\s*$/);
		my ($markerID, $chr, $pos, $marker_type, @sample_types) = split ;
		if ($markerID =~ /ID/){
			@sample_name = @sample_types;
			next ;
		}
		for (my $i=0; $i<@sample_types; $i++){
			my $type = &swap_type($sample_types[$i]);
			push @{$ahtype->{$sample_name[$i]}{$chr}}, $type ;
		}
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

#&stat_marker_info($outfile, \%htype, $asample_name);
sub stat_marker_info()
{
	my ($outfile, $ahtype, $asample_name) = @_ ;
	
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "#Sample\tChr\tTotal_marker\tPresence_marker\tMis_marker\tintegrity\taa_marker\tab_marker\tbb_marker\taa_percent\tab_percent\tbb_percent\ttotal_block\taa_block\tab_block\tbb_block\taa_block_len\tab_block_len\tbb_block_len\n" ;
	for my $sample (@{$asample_name}){
		for my $chr (sort keys %{$ahtype->{$sample}}){
			my @types = @{$ahtype->{$sample}{$chr}};
			my $total_marker = @types ;
			my %hstat = ();
			my %hblock = ();
			my %hblock_len = ();
			my $pre_type = '' ;
			for (my $i=0; $i<@types; $i++){
				$hstat{$types[$i]}++ ;
				next if ($types[$i] eq '--');
				if ($types[$i] ne $pre_type){
					$hblock{$types[$i]}++ ;
					$pre_type = $types[$i] ;
					push @{$hblock_len{$types[$i]}}, 1 ;
				}
				else{
					$hblock_len{$types[$i]}->[-1]++ ;
				}
			}
			my $aa_marker = defined $hstat{aa} ? $hstat{aa} : 0 ;
			my $ab_marker = defined $hstat{ab} ? $hstat{ab} : 0 ;
			my $bb_marker = defined $hstat{bb} ? $hstat{bb} : 0 ;
			my $mis_marker = defined $hstat{'--'} ? $hstat{'--'} : 0 ;
			my $presence_marker = $aa_marker + $ab_marker + $bb_marker ;
			my $integrity = $presence_marker/$total_marker ;
			my $aa_percent = $aa_marker/$presence_marker ;
			my $ab_percent = $ab_marker/$presence_marker ;
			my $bb_percent = $bb_marker/$presence_marker ;
			my $aa_block = defined $hblock{aa} ? $hblock{aa} : 0 ;
			my $ab_block = defined $hblock{ab} ? $hblock{ab} : 0 ;
			my $bb_block = defined $hblock{bb} ? $hblock{bb} : 0 ;
			my $total_block = $aa_block + $ab_block + $bb_block ;
			my $aa_block_len = defined $hblock_len{aa} ? join(',',@{$hblock_len{aa}}) : 0 ;
			my $ab_block_len = defined $hblock_len{ab} ? join(',',@{$hblock_len{ab}}) : 0 ;
			my $bb_block_len = defined $hblock_len{bb} ? join(',',@{$hblock_len{bb}}) : 0 ;
			print OUT "$sample\t$chr\t$total_marker\t$presence_marker\t$mis_marker\t$integrity\t$aa_marker\t$ab_marker\t$bb_marker\t$aa_percent\t$ab_percent\t$bb_percent\t$total_block\t$aa_block\t$ab_block\t$bb_block\t$aa_block_len\t$ab_block_len\t$bb_block_len\n" ;
		}
	}
	close(OUT);

	return ;
}

