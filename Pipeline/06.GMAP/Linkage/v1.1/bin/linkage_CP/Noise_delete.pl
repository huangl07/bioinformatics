#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Storable;
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
my($fout,$floc,$gap,$key,@marker,%marker,%del,%del_num,%rec,$noise_n,$pwd_male,$pwd_female,$map_female,$map_male,$mode,$want_map);
GetOptions(
			"help|?" =>\&USAGE,
			"loc:s"=>\$floc,
			"pwd_male:s"=>\$pwd_male,
			"pwd_female:s"=>\$pwd_female,
			"m_female:s"=>\$map_female,
			"m_male:s"=>\$map_male,
			"o:s"=>\$fout,
			"gap:f"=>\$gap,
			"noise:i"=>\$noise_n,
			"key:s"=>\$key,
			"mode:s"=>\$mode,
			"map:s"=>\$want_map,
		) or &USAGE;
&USAGE unless ($floc and $pwd_male and $pwd_female and $map_female and $map_male and $fout);

$gap||=0.1;
$noise_n||=1;
my $gap_x=0.01;
$mode=0 if(!defined $mode);
my($map,$pwd,$fkey);
if ($mode==0) {
	$map=$map_male;
	$pwd=$pwd_male;
	$fkey=$key."male";
}else{
	$map=$map_female;
	$pwd=$pwd_female;
	$fkey=$key."female";
}

$floc=AbsolutePath("file",$floc) if($floc);
my $floc_base=basename($floc) if($floc);
$map=AbsolutePath("file",$map);
my $map_base=basename($map);
$pwd=AbsolutePath("file",$pwd);
my $pwd_base=basename($pwd);
mkdir "$fout" unless(-d $fout);
$fout=AbsolutePath("dir",$fout);
mkdir "$fout/data" if(!-d "$fout/data");
my $step=0;

if ($step==0) {
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


	open ON,"$pwd";
	while (<ON>) {
		chomp;
		next if(/^\s+$/ or /^$/ or $_!~/Marker/i);
		my($marker1,$marker2,$rec,$lod)=split /\s+/,$_;
		$rec{$marker1}{$marker1}{'r'}=0;
		$rec{$marker2}{$marker2}{'r'}=0;
		$rec{$marker1}{$marker2}{'r'}=$rec;
		$rec{$marker2}{$marker1}{'r'}=$rec;
		$rec{$marker1}{$marker1}{'lod'}=0;
		$rec{$marker2}{$marker2}{'lod'}=0;
		$rec{$marker1}{$marker2}{'lod'}=$lod;
		$rec{$marker2}{$marker1}{'lod'}=$lod;
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


	my %del_marker;
	foreach my $m1 (keys %del) {
		foreach my $m2 (keys %{$del{$m1}}) {
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
		print OUT1 "$m\t$marker{$m}\n";
		for (my $i=$j;$i<@final_marker;$i++) {
			print OUT2 "$m\t$final_marker[$i]\t$rec{$m}{$final_marker[$i]}{'r'}\t$rec{$m}{$final_marker[$i]}{'lod'}\n";
		}
		$j+=1;
	}
	close OUT1;
	close OUT2;
	my $draw_heatmap = "perl $Bin/drawFlow/drawHeatmap.pl -p $fout/data/correct\_$pwd_base -i $fout/data/correct\_$map_base -k  $fkey -d $fout/map" ;
	`$draw_heatmap`;


	if (glob "$fout/data/*.array") {
		my($array_data)=glob "$fout/data/*.array";
		my $array=retrieve("$array_data");
		unlink "$array_data";
		foreach (@$array){
			$del_marker{$_}=1;
		}
		my $marker_num=`less -S $floc \| grep \'Marker\' \| wc -l`;
		$marker_num-=scalar(keys %del_marker);
		open ON,"$floc" || die "please input .loc";
		open OUT,">$fout/data/correct\_$floc_base";

		while (<ON>) {
			if($_!~/Marker/i && $_!~/nloc/i){
				print OUT "$_";
			}elsif (/nloc/i) {
				print OUT "nloc = $marker_num\n"
			}else{
				my($marker)=split /\s+/,$_;
				if(grep{$marker eq $_}keys %del_marker){
					next;
				}else{
					print OUT $_ ;
				}
			}
		}
		close ON;
		close OUT;
		exit(1);
	}else{
		my $array=[keys %del_marker];
		store  $array,"$fout/data/$fkey.array";
		my $command="perl $0 -loc $floc  -pwd_male $pwd_male -pwd_female $pwd_female -m_female $map_female -m_male $map_male -o $fout -gap $gap -noise $noise_n  -mode 1 -key $key";
		`$command`;
	}
}
#
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
	杂点清除，然后再排图
	Options:
			"help|?" =>\&USAGE,
			"loc:s"=>\$floc,       must  all the loc
			"pwd_male:s"=>\$pwd_male,
			"pwd_female:s"=>\$pwd_female,
			"m_female:s"=>\$map_female,
			"m_male:s"=>\$map_male,
			"o:s"=>\$fout,
			"gap:f"=>\$gap,
			"noise:i"=>\$noise_n,
			"key:s"=>\$key,
			"mode:s"=>\$mode,
USAGE
	print $usage;
	exit;
}