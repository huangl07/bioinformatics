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
GetOptions(\%opts,"i=s","od=s","p=s","w=s","s=s","ad=s","bd=s","h" );

#&help()if(defined $opts{h});
if(!defined($opts{i}) || !defined($opts{od}) || !defined($opts{p}) || !defined($opts{w}) || !defined($opts{s}) || !defined($opts{ad}) || !defined($opts{bd}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		This script is write to imputation marker genotype.
		Version: $ver
		v1.0:	Function realization.
		v1.1:	Memory optimization.

	Usage:
		-i           genotype file               <infile>     must be given
		-od          result file dir             <outdir>     must be given
		-p           outfile prefix              <string>     must be given
		-w           window length               <int>        must be given
		-s           step length                 <int>        must be given
		-ad          a depth for judge as aa     <int>        must be given
		-bd          b depth for judge as bb     <int>        must be given
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
my $outdir = $opts{od} ;
mkdir $outdir if (!-d $outdir) ;
$outdir = &ABSOLUTE_DIR($outdir);
my $prefix = $opts{p} ;
my $win_len = $opts{w} ;
my $step_len = $opts{s} ;
my $a_depth = $opts{ad} ;
my $b_depth = $opts{bd} ;

## reading genotype file and store
&show_log("reading genotype file.");
my %hmarker = ();
my %hmarker_info = ();
my %hmarker_num = ();
my $asample_name = &get_marker_info($marker_file, \%hmarker, \%hmarker_info, \%hmarker_num) ;
&show_log("reading genotype file done.");

## slide window and get type and imputation
&show_log("get window type.");
my %himp_marker = ();
&slide_window(\%hmarker, \%hmarker_num, \%himp_marker);
&show_log("get window type done.");

## imputation for missing marker
#&show_log("imputation missing marker.");
#my %himp_marker = ();
#&imputation(\%hwindow, \%hmarker_num, \%himp_marker);
#&show_log("imputation missing marker done.");

## out put
&show_log("out put marker.");
&output_marker(\%himp_marker, \%hmarker_info, $asample_name, $outdir, $prefix);
&show_log("out put marker done.");


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

#my $asample_name = &get_marker_info($marker_file, \%hmarker, \%hmarker_info, \%hmarker_num) ;
sub get_marker_info()
{
	my ($marker_file, $ahmarker, $ahmarker_info, $ahmarker_num) = @_ ;

	open (IN, $marker_file) || die "Can't creat $marker_file, $!\n" ;
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

#&slide_window(\%hmarker, \%hmarker_num, \%himp_marker);
sub slide_window()
{
	my ($ahmarker, $ahmarker_num, $ahimp_marker) = @_ ;

	for my $sample (keys %{$ahmarker}){
		for my $chr (keys %{$ahmarker->{$sample}}){
			my $amarker = &get_window_marker($ahmarker->{$sample}->{$chr});
			my $aimp_marker = &imput_chr_marker($amarker, $ahmarker_num->{$chr});
			$ahimp_marker->{$sample}{$chr} = $aimp_marker ;
		}
	}

	return ;
}

#my $amarker = &get_window_marker();
sub get_window_marker()
{
	my ($atypes) = @_ ;

	my $last_flag = 0 ;
	my @marker = ();
	for (my $i=0; $i<@{$atypes}; $i++){
		next if ($atypes->[$i] eq '--');
		my $num = 0 ;
		my %htypes = ();
		for (my $j=$i; $j<@{$atypes}; $j++){
			next if ($atypes->[$j] eq '--');
			$htypes{$atypes->[$j]}++ ;
			$num++ ;
			$last_flag = 1 if ($j == $#$atypes) ;
			last if ($num % $win_len == 0);
		}
		last if ($num % $win_len != 0);
		my $win_type = &judge_type(\%htypes);
		$marker[$i] = $win_type ;
		last if ($last_flag == 1);
	}

	return(\@marker);
}

#my $win_type = &judge_type(\%htypes);
sub judge_type()
{
	my ($ahtype) = @_ ;

	my $aa_num = defined $ahtype->{'aa'} ? $ahtype->{'aa'} : 0 ;
	my $ab_num = defined $ahtype->{'ab'} ? $ahtype->{'ab'} : 0 ;
	my $bb_num = defined $ahtype->{'bb'} ? $ahtype->{'bb'} : 0 ;

	my $a_num = $aa_num*2 + $ab_num ;
	my $b_num = $bb_num*2 + $ab_num ;

	if ($a_num >= 2*$a_depth){
		return ('aa');
	}
	elsif ($b_num >= 2*$b_depth){
		return ('bb');
	}
	else {
		return ('ab');
	}
}

#&imputation(\%hwindow, \%hmarker_num, \%himp_marker);
sub imputation()
{
	my ($ahwindow, $ahmarker_num, $ahimp_marker) = @_ ;

	for my $sample (keys %{$ahwindow}){
		for my $chr (keys %{$ahwindow->{$sample}}){
			my $aimp_marker = &imput_chr_marker($ahwindow->{$sample}->{$chr}, $ahmarker_num->{$chr});
			$ahimp_marker->{$sample}{$chr} = $aimp_marker ;
		}
	}

	return ;
}

#my $aimp_marker = &imput_chr_marker($ahwindow->{$sample}->{$chr}, $marker_num);
sub imput_chr_marker()
{
	my ($awindow_marker, $marker_num) = @_ ;
	my @imp_marker = ();

	my ($pre_marker, $pre_index) = &find_next_marker($awindow_marker, 0);
	if ($pre_index == -1){
		for (my $i=0; $i<$marker_num; $i++){
			$imp_marker[$i] = '--' ;
		}
		return(\@imp_marker);
	}
	for (my $i=0; $i<$pre_index; $i++){
		$imp_marker[$i] = $pre_marker ;
	}
	my ($next_marker, $next_index) = &find_next_marker($awindow_marker, $pre_index+1);
	for (my $i=$pre_index; $i<$marker_num; $i++){
		if (!defined $awindow_marker->[$i]){
			if ($pre_marker eq $next_marker){
				$imp_marker[$i] = $pre_marker ;
			}
			else{
				$imp_marker[$i] = '--' ;
			}
		}
		else{
			$imp_marker[$i] = $awindow_marker->[$i] ;
			$pre_marker = $awindow_marker->[$i] ;
			$pre_index = $i ;
			($next_marker, $next_index) = &find_next_marker($awindow_marker, $pre_index+1);
			last if ($next_marker eq '-1') ;
		}
	}
	for (my $i=$pre_index; $i<$marker_num; $i++){
		$imp_marker[$i] = $pre_marker ;
	}

	return(\@imp_marker);
}

#my ($marker, $index) = &find_next_marker($awindow_marker, $start_index);;
sub find_next_marker()
{
	my ($awindow_marker, $start_index) = @_ ;

	for (my $i=$start_index; $i<@{$awindow_marker}; $i++){
		if (defined $awindow_marker->[$i]){
			return ($awindow_marker->[$i], $i);
		}
	}

	return('-1', '-1');
}

#&output_marker(\%himp_marker, \%hmarker_info, $asample_name, $outdir, $prefix);
sub output_marker()
{
	my ($ahimp_marker, $ahmarker_info, $asample_name, $outdir, $prefix) = @_ ;

	my $outfile = "$outdir/$prefix.imp" ;
	open (OUT, ">$outfile") || die "Can't creat $outfile, $!\n" ;
	print OUT "MarkerID\tChr\tpos\ttype\tP\tM\t", join("\t", @{$asample_name}), "\n" ;
	for my $chr (sort keys %{$ahmarker_info}){
		for (my $i=0; $i<@{$ahmarker_info->{$chr}}; $i++){
			print OUT $ahmarker_info->{$chr}[$i][0], "\t$chr\t", $ahmarker_info->{$chr}[$i][1],"\taaxbb\taa\tbb" ;
			for my $sample (@{$asample_name}){
				print OUT "\t", $ahimp_marker->{$sample}{$chr}->[$i] ;
			}
			print OUT "\n" ;
		}
	}
	close(OUT);

	return ;
}


