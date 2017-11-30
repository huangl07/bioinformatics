#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use POSIX qw(strftime);
my $BEGIN_TIME=time();
my $cur_time=strftime("%Y-%m-%d %A %H:%M:%S", localtime());
my $version="1.0";
$Script=~s/\.pl//g;
my($fIn,$pwd,$fout,$map,$floc,$gap,$fkey,@marker,%marker,%del,%del_num,%rec,$noise_n,$popt,$want_map);
GetOptions(
			"help|?" =>\&USAGE,
			"i:s"=> \$fIn,
			"pwd:s"=>\$pwd,
			"m:s"=>\$map,
			"loc:s"=>\$floc,
			"o:s"=>\$fout,
			"gap:s"=>\$gap,
			"noise:i"=>\$noise_n,
			"key:s"=>\$fkey,
			"popt:s"=>\$popt,
			"map:s"=>\$want_map,
		) or &USAGE;
if($popt=~/BC1/i){
	&USAGE unless ($pwd);
}else{
	&USAGE unless ($fIn);
}
&USAGE unless ($map and $fout and $floc and $popt);
$fkey||='result';
$gap||=0.1;
$noise_n||=1;
my $gap_x=0.01;

$floc=AbsolutePath("file",$floc);
my $floc_base=basename($floc);
$map=AbsolutePath("file",$map);
my $map_base=basename($map);
$fIn=AbsolutePath("file",$fIn) if($fIn);
my $fIn_base=basename($fIn) if($fIn);
$pwd=AbsolutePath("file",$pwd) if($pwd);
my $pwd_base=basename($pwd) if($pwd);
mkdir "$fout" unless(-d $fout);
$fout=AbsolutePath("dir",$fout);
mkdir "$fout/data" if(!-d "$fout/data");
open ON,"$map";
my $map_head;
while (<ON>) {
	if ($_!~/Marker/i || /^\s+$/){
		$map_head.=$_;
		next;
	}
	chomp;
	my($marker,$posi)=split /\s+/,$_;
	$marker{$marker}=$posi;
	push @marker,$marker;
}
close ON;
if ($popt=~/bc\d+/i || $popt=~/DH/i) {
	open ON,"$pwd";
	while (<ON>) {
		chomp;
		next if(/^\s+$/ or /^$/ or $_!~/Marker/i);
		my($marker1,$marker2,$rec,$mlod)=split /\s+/,$_;
		$rec{$marker1}{$marker1}{'r'}=0;
		$rec{$marker2}{$marker2}{'r'}=0;
		$rec{$marker1}{$marker1}{'lod'}=0;
		$rec{$marker2}{$marker2}{'lod'}=0;
		$rec{$marker1}{$marker2}{'r'}=$rec;
		$rec{$marker2}{$marker1}{'r'}=$rec;
		$rec{$marker1}{$marker2}{'lod'}=$mlod;
		$rec{$marker2}{$marker1}{'lod'}=$mlod;
	}
	close ON;
	my $n=0;
	foreach my $marker (@marker) {
		my @marker_cycle=@marker;
		next if(!exists $rec{$marker});
		my $i=$n+1;
		while (1) {
			last if ($i>=@marker-2);
			if($marker_cycle[$i] eq 'X'){
				$i++;next;
			}
			if(!exists $rec{$marker}{$marker_cycle[$i]}){
				$i++;next;
			}
			if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_right(\@marker_cycle,$i-1,$marker)}{'r'})<$gap) {
				$i++;
			}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_right(\@marker_cycle,$i-1,$marker)}{'r'})>=$gap){
				if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i+1]}{'r'})>=0.95*$gap) {
					$del_num{$marker}++;$del_num{$marker_cycle[$i]}++;
					$del{$marker}{$marker_cycle[$i]}=1;$del{$marker_cycle[$i]}{$marker}=1;
					$marker_cycle[$i]='X';
					$i++;
				}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i+1]}{'r'})<$gap_x) {
					my $j=$i+1;
					while (1) {
						if ($j==@marker_cycle-1){
							last;
						}
						if (abs($rec{$marker}{$marker_cycle[$j]}{'r'}-$rec{$marker}{$marker_cycle[$j-1]}{'r'})<$gap_x || $rec{$marker}{$marker_cycle[$j]}{'r'}>=$rec{$marker}{$marker_cycle[$i]}{'r'}){
							$j++;
						}else{
							last;
						}
						
					}
					if ($j-$i<=$noise_n && $j!=@marker_cycle-1) {
						foreach  (@marker_cycle[$i..$j]) {
							$del_num{$marker}++;$del_num{$_}++;
							$del{$marker}{$_}=1;$del{$_}{$marker}=1;
						}
						@marker_cycle[$i..$j]=split //,('X' x ($j-$i+1));
					}
					$i=$j+1;
				}else{
					$i++;
				}
			}
		}

		$i=$n-1;
		while (1) {
			last if ($i<=1);
			if($marker_cycle[$i] eq 'X'){
				$i--;next;
			}
			if(!exists $rec{$marker}{$marker_cycle[$i]}){
				$i--;next;
			}
			if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_left(\@marker_cycle,$i+1,$marker)}{'r'})<$gap) {
				$i--;
			}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_left(\@marker_cycle,$i+1,$marker)}{'r'})>=$gap){
				if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i-1]}{'r'})>=0.95*$gap) {
					$del_num{$marker}++;$del_num{$marker_cycle[$i]}++;
					$del{$marker}{$marker_cycle[$i]}=1;$del{$marker_cycle[$i]}{$marker}=1;
					$marker_cycle[$i]='X';
					$i--;
				}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i-1]}{'r'})<$gap_x) {
					my $j=$i-1;
					while (1) {
						if ($j==0){
							last;
						}
						if (abs($rec{$marker}{$marker_cycle[$j]}{'r'}-$rec{$marker}{$marker_cycle[$j+1]}{'r'})<$gap_x || $rec{$marker}{$marker_cycle[$j]}{'r'}>=$rec{$marker}{$marker_cycle[$i]}{'r'}){
							$j--;
						}else{
							last;
						}
						
					}
					if ($i-$j<=$noise_n && $j!=0) {
						foreach  (@marker_cycle[$j..$i]) {
							$del_num{$marker}++;$del_num{$_}++;
							$del{$marker}{$_}=1;$del{$_}{$marker}=1;
						}
						@marker_cycle[$j..$i]=split //,('X' x ($i-$j+1));
					}
					$i=$j-1;
				}else{
					$i--;
				}
			}
		}
		$n+=1;
	}
	if(keys %del==0){
		warn "their are not some noise!";
		exit(-1);
	}
	
	my %del_marker;
	foreach my $m1 (keys %del) {
		foreach my $m2 (keys %{$del{$m1}}) {
		#	next if($del_num{$m1}<=0 || $del_num{$m2}<=0);
		#	next if($del_num{$m2}<=0);
			my $delete=$del_num{$m1}>$del_num{$m2}?$m1:$m2;
			$del_marker{$delete}=1;
		}
	}



	my @final_marker;
	foreach my $x (@marker) {
		next if(grep{$x eq $_}keys %del_marker);
		push @final_marker,$x;
	}
	open OUT1,">$fout/data/correct\_$map_base";
	open OUT2,">$fout/data/correct\_$pwd_base";
	print OUT1 "$map_head";
	my $j=1;
	foreach my $m (@final_marker) {
		print OUT1 "$m\t$marker{$m}\t\n";
		for (my $i=$j;$i<@final_marker;$i++) {
			print OUT2 "$m\t$final_marker[$i]\t$rec{$m}{$final_marker[$i]}{'r'}\t$rec{$m}{$final_marker[$i]}{'lod'}\n";
		}
		$j+=1;
	}
	close OUT1;
	close OUT2;
	my $draw_heatmap = "perl $Bin/draw/draw_heatmap.pl -p $fout/data/correct\_$pwd_base -i $fout/data/correct\_$map_base -k  $fkey -d $fout/map" ;
	`$draw_heatmap`;
	open ON,"$floc";
	open OUT,">$fout/data/correct\_$floc_base";
	while (<ON>) {
		if($_!~/Marker/i){
			print OUT $_;
			next;
		}
		my($marker)=split /\s+/,$_;
		print OUT $_ if(grep{$marker eq $_}@final_marker);
	}
	close ON;
	close OUT;
}else{
	my $n=-2;
	open ON,"$fIn" || die "$!\n";
	#open OUT,">./var.csv";
	my $rec_head;
	while (<ON>) {
		$rec_head.=$_ if($n<0);
		chomp;
		if($n==-2 && /^=+$/){
			$n+=1;
		}elsif($n==-1 && /^=+$/){
			$n+=1; 
		}elsif ($n==0 && /marker/i) {
			$n+=1; 
		}elsif ($n>0) {
			my($marker,@rec)=split /\s+/,$_;
			for (my $i=0;$i<@marker;$i++) {
				$rec{$marker}{$marker[$i]}{'r'}=$rec[$i];
				$rec{$marker[$i]}{$marker}{'r'}=$rec[$i];
			}
		}
	}
	close ON;
	$n=0;
	foreach my $marker (@marker) {
		my @marker_cycle=@marker;
		next if(!exists $rec{$marker});
		my $i=$n+1;
		while (1) {
			last if ($i>=@marker-2);
			if($marker_cycle[$i] eq 'X'){
				$i++;next;
			}
			if(!exists $rec{$marker}{$marker_cycle[$i]}){
				$i++;next;
			}
			if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_right(\@marker_cycle,$i-1,$marker)}{'r'})<$gap) {
				$i++;
			}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_right(\@marker_cycle,$i-1,$marker)}{'r'})>=$gap){
				if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i+1]}{'r'})>=0.95*$gap) {
					$del_num{$marker}++;$del_num{$marker_cycle[$i]}++;
					$del{$marker}{$marker_cycle[$i]}=1;$del{$marker_cycle[$i]}{$marker}=1;
					$marker_cycle[$i]='X';
					$i++;
				}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i+1]}{'r'})<$gap_x) {
					my $j=$i+1;
					while (1) {
						if ($j==@marker_cycle-1){
							last;
						}
						if (abs($rec{$marker}{$marker_cycle[$j]}{'r'}-$rec{$marker}{$marker_cycle[$j-1]}{'r'})<$gap_x || $rec{$marker}{$marker_cycle[$j]}{'r'}>=$rec{$marker}{$marker_cycle[$i]}{'r'}){
							$j++;
						}else{
							last;
						}
						
					}
					if ($j-$i<=$noise_n && $j!=@marker_cycle-1) {
						foreach  (@marker_cycle[$i..$j]) {
							$del_num{$marker}++;$del_num{$_}++;
							$del{$marker}{$_}=1;$del{$_}{$marker}=1;
						}
						@marker_cycle[$i..$j]=split //,('X' x ($j-$i+1));
					}
					$i=$j+1;
				}else{
					$i++;
				}
			}
		}

		$i=$n-1;
		while (1) {
			last if ($i<=1);
			if($marker_cycle[$i] eq 'X'){
				$i--;next;
			}
			if(!exists $rec{$marker}{$marker_cycle[$i]}){
				$i--;next;
			}
			if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_left(\@marker_cycle,$i+1,$marker)}{'r'})<$gap) {
				$i--;
			}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{&marker_r_left(\@marker_cycle,$i+1,$marker)}{'r'})>=$gap){
				if (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i-1]}{'r'})>=0.95*$gap) {
					$del_num{$marker}++;$del_num{$marker_cycle[$i]}++;
					$del{$marker}{$marker_cycle[$i]}=1;$del{$marker_cycle[$i]}{$marker}=1;
					$marker_cycle[$i]='X';
					$i--;
				}elsif (abs($rec{$marker}{$marker_cycle[$i]}{'r'}-$rec{$marker}{$marker_cycle[$i-1]}{'r'})<$gap_x) {
					my $j=$i-1;
					while (1) {
						if ($j==0){
							last;
						}
						if (abs($rec{$marker}{$marker_cycle[$j]}{'r'}-$rec{$marker}{$marker_cycle[$j+1]}{'r'})<$gap_x || $rec{$marker}{$marker_cycle[$j]}{'r'}>=$rec{$marker}{$marker_cycle[$i]}{'r'}){
							$j--;
						}else{
							last;
						}
						
					}
					if ($i-$j<=$noise_n && $j!=0) {
						foreach  (@marker_cycle[$j..$i]) {
							$del_num{$marker}++;$del_num{$_}++;
							$del{$marker}{$_}=1;$del{$_}{$marker}=1;
						}
						@marker_cycle[$j..$i]=split //,('X' x ($i-$j+1));
					}
					$i=$j-1;
				}else{
					$i--;
				}
			}
		}
		$n+=1;
	}


	if(keys %del==0){
		warn "their are not some noise!";
		exit(-1);
	}
	my @del_marker;
	foreach my $m1 (keys %del) {
		foreach my $m2 (keys %{$del{$m1}}) {
		#	next if($del_num{$m1}<=0 || $del_num{$m2}<=0);
			my $delete=$del_num{$m1}>$del_num{$m2}?$m1:$m2;
			push @del_marker,$delete;
		}
	}
	my @final_marker;
	foreach my $x (@marker) {
		next if(grep{$x eq $_}@del_marker);
		push @final_marker,$x;
	}
	
	open OUT1,">$fout/data/correct\_$map_base";
	open OUT2,">$fout/data/correct\_$fIn_base";
	print OUT1 "$map_head";
	print OUT2 "$rec_head";
	print OUT2 "\t",(join "\t",@final_marker),"\n";
	foreach my $m1 (@final_marker) {
		print OUT1 "$m1\t$marker{$m1}\n";
		print OUT2 "$m1\t";
		my @line;
		foreach my $m2 (@final_marker) {
			push @line,$rec{$m1}{$m2};
		}
		print OUT2 (join "\t",@line),"\n";
		
	}
	close OUT1;
	close OUT2;
	my $draw_heatmap = " perl $Bin/draw/draw_MST_heatmap.pl -i $fout/data/correct\_$fIn_base -k  $fkey -d $fout/map";
	`$draw_heatmap`;
	open ON,"$floc";
	open OUT,">$fout/data/correct\_$floc_base";
	while (<ON>) {
		if($_!~/Marker/i && $_!~/nloc=/i){
			print OUT $_;
		}elsif (/nloc=/i) {
			my $nloc=@final_marker;
			print OUT "nloc = $nloc\n"
		}else{
			my($marker)=split /\s+/,$_;
			print OUT $_ if(grep{$marker eq $_}@final_marker);
		}
	}
	close ON;
	close OUT;
}

if(!defined $want_map){
	print "Do you want mapping? (yes or no)\n";
	$want_map=<STDIN>;
	#while ($want_map=~/^\s+$/ || $want_map eq '' || $want_map!~/yes/ || $want_map!~/no/) {
	while($want_map!~/\byes\b/ && $want_map!~/\bno\b/){
		$want_map=<STDIN>;
	}
	if ($want_map=~/yes/) {
		my($chr)=$floc_base=~/^([^\.]+)\..+\./;
		my $command="perl $Bin/sgsMapSmoothCycle.pl -i $fout/data/correct\_$floc_base -k $chr -d $fout/correct_map -mode 0 -n 1 ";
		print "mapping is start...\n";
		`$command`;
	}
}else{
	my($chr)=$floc_base=~/^([^\.]+)\..+\./;
	my $command="perl $Bin/sgsMapSmoothCycle.pl -i $fout/data/correct\_$floc_base -k $chr -d $fout/correct_map -mode 0 -n 1 ";
	print "mapping is start...\n";
	`$command`;
}



########################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\t",$cur_time,"\n";
########################################################################################


# -----------------------------------------------------------sub function---------------------------------------------------#
sub marker_r_right{
	my($marker_array,$cur_i,$marker)=@_;
	while ($marker_array->[$cur_i] eq 'X') {
		$cur_i--;
		last if($marker_array->[$cur_i] eq $marker);
	}
	return $marker_array->[$cur_i];

}
sub marker_r_left{
	my($marker_array,$cur_i,$marker)=@_;
	while ($marker_array->[$cur_i] eq 'X') {
		$cur_i++;
		last if($marker_array->[$cur_i] eq $marker);
	}
	return $marker_array->[$cur_i];

}

sub AbsolutePath{
        my ($type,$input) = @_;
        my $return;
        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chop $return;
                $return .="/".$file;
                chdir($pwd);
        }
        return $return;
}


# -----------------------------------------------------------sub USAGE---------------------------------------------------#
sub USAGE {
    my $usage=<<"USAGE";
    Program:$Script
    Version:$version    $cur_time
    Contact:zhaogf <zhaogf\@biomarker.com.cn>
	Options:
		"help|?"
		"i:s"		.mst.o (Ri f2)
		"pwd:s"		.pwd (BC1 DH)
		"m:s"		.map
		"loc:s"		.loc
		"o:s"        outdir
		"gap:s"      if great than gap, then the point is noise  default 0.1 
		"noise:i"	 number of noise default 1
		"key:s"
		"popt:s"   
USAGE
	print $usage;
	exit;
}