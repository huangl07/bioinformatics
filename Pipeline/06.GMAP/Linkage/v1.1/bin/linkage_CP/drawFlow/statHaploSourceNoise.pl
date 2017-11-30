#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$outdir,$sus_ratio,$file_kind);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
#				"id:s"=>\$indir_haploSource,
#				"k:s"=>\$fKey,
				"od:s"=>\$outdir,
				"r:f"=>\$sus_ratio,
				"kind:s"=>\$file_kind,
				) or &USAGE;
&USAGE unless ($fIn );

$sus_ratio||=0.06;
$file_kind||="sexAver";
#-------------------------------------------------------------------
#Global parameter settings
#-------------------------------------------------------------------
$outdir||="StatHaploSourceNoise";
mkdir $outdir unless (-d $outdir) ;
$outdir=AbsolutePath("dir",$outdir);

#输入可能是文件，也可能是带有文件的目录
my @file=();
if (-f $fIn) {
	$fIn=AbsolutePath("file",$fIn);
	push @file, $fIn;
}elsif (-d $fIn){
	my $indir=AbsolutePath("dir",$fIn);
	if ($file_kind eq "sexAver") {
		@file=glob("$indir/*.sexAver.phase");
	}elsif($file_kind eq "all"){
		@file=glob("$indir/*.phase");
	}

}


#######  Get xxx.phase file ##############

foreach my $fhaploSource (@file) {
	my (%haplo,@marker,%indiToCheck,%depth,%genotype,%allel,%depthStat,%cross_type,%G,%linkagePhase)=();
	my %suspicious;
	my %sus_marker;
	my %order=();
	my %order_marker=();
	my $order=0;
	######## Get Data
	open (H,"$fhaploSource") or die $!;
	while (<H>) {
		chomp;
		next if (/^$/) ;
		my ($marker,$crossType,$linkPhase,@haploSource)=split;
		for (my $i=0;$i<@haploSource ;$i++) {
			my @unit=split //,$haploSource[$i];
			push @{$haplo{$i}},[@unit];
		}
		push @marker,$marker;
		$order_marker{$order}=$marker;
		$order{$marker}=$order;
		$order++;
	}

	close (H) ;
	my $marker_sum=@marker;
	my $off_sum=scalar keys %haplo;

	####统计marker的可疑子代数
	foreach my $indi (sort {$a <=> $b} keys %haplo) {
		for (my $chain_num=0;$chain_num<2;$chain_num++) {
			for (my $i=0;$i<@{$haplo{$indi}} ;$i++) {
				next if($haplo{$indi}->[$i]->[$chain_num] ne "1" && $haplo{$indi}->[$i]->[$chain_num] ne "2");
				my $before=$haplo{$indi}->[$i]->[$chain_num];
				my $after=$haplo{$indi}->[$i]->[$chain_num];
#				if ($haplo{$indi}->[$i]->[$chain_num] eq "1") {
#					$before="2";
#					$after="2";
#				}else{
#					$before="1";
#					$after="1";
#				}
				my $sign='0';
				for (my $j=$i-1;$j>=0;$j--) {
					if ($haplo{$indi}->[$j]->[$chain_num] eq "1" || $haplo{$indi}->[$j]->[$chain_num] eq "2") {
						if ($i==@{$haplo{$indi}}-1) {
							if ($sign eq "0") {
								$sign=$haplo{$indi}->[$j]->[$chain_num];
							}elsif ($sign eq "2" || $sign eq "1") {
								if ($sign eq $haplo{$indi}->[$j]->[$chain_num]) {
									$before=$haplo{$indi}->[$j]->[$chain_num];
									$after=$haplo{$indi}->[$j]->[$chain_num];
								}
								last;
							}
						}else{
							$before=$haplo{$indi}->[$j]->[$chain_num];
							last;
						}
					}
				}
				$sign='0';
				for (my $j=$i+1;$j<@{$haplo{$indi}};$j++) {
					if ($haplo{$indi}->[$j]->[$chain_num] eq "1" || $haplo{$indi}->[$j]->[$chain_num] eq "2") {
						if ($i==0) {
							if ($sign eq "0") {
								$sign=$haplo{$indi}->[$j]->[$chain_num];
							}elsif ($sign eq "2" || $sign eq "1") {
								if ($sign eq $haplo{$indi}->[$j]->[$chain_num]) {
									$before=$haplo{$indi}->[$j]->[$chain_num];
									$after=$haplo{$indi}->[$j]->[$chain_num];
#									print $indi,"\t",$haplo{$indi}->[$i]->[$chain_num],"\t",$after,"\n";
								}
								last;
							}
						}else{
							$after=$haplo{$indi}->[$j]->[$chain_num];
							last;
						}
					}
				}
				if ($before eq $after && $haplo{$indi}->[$i]->[$chain_num] ne $before ) {
					$suspicious{$marker[$i]}{$indi}=1;
				}
			}
		}
	}

	foreach my $marker (keys %suspicious) {
		my $sus_off_sum=scalar keys %{$suspicious{$marker}};
		if ($sus_off_sum>$off_sum*$sus_ratio) {
			$sus_marker{$marker}=$sus_off_sum;
		}
	}

	#-------------------------------------------------------------------
	# Print
	#-------------------------------------------------------------------
	my $noise_total_sum=0;
	foreach my $marker (keys %suspicious) {
		my $sus_off_sum=scalar keys %{$suspicious{$marker}};
		$noise_total_sum+=$sus_off_sum;
	}
	my $gene_sum=$off_sum*$marker_sum;
	my $noise_ratio=$noise_total_sum/$gene_sum;

	my $fhaploSource_name=basename($fhaploSource);
	open (SUS,">","$outdir/$fhaploSource_name.noiseStat") or die $!;
	print SUS "marker sum : $marker_sum\n";
	print SUS "ind sum : $off_sum\n";
	print SUS "gene sum : $gene_sum\n";
	print SUS "noise sum : $noise_total_sum\n";
	print SUS "noise ratio : ",sprintf("%.3f",$noise_ratio*100),"%\n";
	print SUS "\n";
	print SUS "ratio : $sus_ratio\n";
	print SUS "noise threshold :",$off_sum*$sus_ratio,"\n";
#	print SUS "\n";
	my $sus_marker_sum=keys %sus_marker;
	print SUS "Suspicious marker sum : $sus_marker_sum\n\n";
	print SUS "Marker\tNoise_Sum\tOrder_In_Map\n";
	foreach my $marker (sort {keys %{$suspicious{$b}} <=> keys %{$suspicious{$a}} } keys %sus_marker ) {
		my $sus_sum=scalar keys %{$suspicious{$marker}};
		print SUS $marker,"\t",$sus_sum,"\t",$order{$marker}+1,"\n";
	}
	print SUS "\n";
	print SUS "****************************************************************\n";
	print SUS "Marker\tNoise_Sum\n";
	foreach my $marker (@marker) {
		my $sus_off_sum=0;
		if (exists $suspicious{$marker}) {
			$sus_off_sum=scalar keys %{$suspicious{$marker}};
		}
		print SUS $marker,"\t",$sus_off_sum,"\n";
	}
	close (SUS) ;
}




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub AbsolutePath
{		#获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
	$/="\n";

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
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Ma chouxian <macx\@biomarker.com.cn> 
Discription: 
	stat noise point of haploSource pictrue.
	creat dir : key.StatHaploSourceNoise/.

Usage:
  Options:
  -i		<file>	inputfile or dir of inputfile(xxx.phase)                                   forced
  -od		<str>	dir of output file,                                                        default key.StatHaploSourceNoise
  -r		<float>	ratio threshold between noise sum and total indi sum,                      default 0.06
  -kind		<str>	"sexAver" indicate only deal with file with name of "*sexAver.phase",
                    "all" indicate deal with file of "*.phase"                                 default "sexAver"
  -h		Help

USAGE
	print $usage;
	exit;
}