#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Storable;
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($dIn,$fKey,$dOut,$mapFunction,$gapThreshold,$rThreshold,$lodThreshold,$poorFitThreshold);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$dIn,
				"k:s"=>\$fKey,
				"d:s"=>\$dOut, 
				
				"gap:s"=>\$gapThreshold, 
				"mapFunctiion:s"=>\$mapFunction, 
				"r:s"=>\$rThreshold, 
				"lod:s"=>\$lodThreshold, 	

				) or &USAGE;
&USAGE unless ($dIn and $fKey);

# ------------------------------------------------------------------
# global parameter's settings 
# ------------------------------------------------------------------

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;

$mapFunction||= "Kosambi";
$gapThreshold||=5;
$rThreshold||=0.4;
$lodThreshold||=1;

# ------------------------------------------------------------------
# global value
# ------------------------------------------------------------------

my (%Map,%loc,%pwd,%EvaluationStat) = ();

my %mapFunction=(
	"Kosambi"=>\&Kosambi,
	"Haldane"=>\&Haldane,
	);
my %inverseMapFunction=(
	"Kosambi"=>\&inverseKosambi,
	"Haldane"=>\&inverseHaldane,
	);
# ------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------
&loadMap(\%Map,$dIn);
&loadLOC(\%loc,$dIn);
&loadPWD(\%pwd,$dIn);

#print Dumper %Map;die;
#print Dumper %loc;die;
#print Dumper %pwd;die;

# ------------------------------------------------------------------
# statistic
# ------------------------------------------------------------------

# divide the evaluation into following parts:

# basic statistic of maps (sexAver female male ),include the number of loci,total map distance,average map distance,ratio of gap less than five centiMorgans
# final chisquare,RMSD,Average map error to evaluate the entire map quality
# details of marker : error genotype,N.N.fit,mean chisquare contribution,jump 

print "statistic ...\n";
open (OUT,">$dOut/$fKey.map.stat") or die $!;
print OUT ";line starts with '>' : linkage_group_ID\tsex\tmarker_num\ttotal_CM\taver_CM\tgap_lt$gapThreshold\tchisquare\tdf\tmean_chisquare\tRMSD\tAFE\n";
print OUT ";line starts with LG : map_order\tmarkerID\ttype\tmap_position\tN.N.Fit_rUnit\tN.N.Fit_cmUnit\tr_synteny\tjump\n";
foreach my $lg (sort { $a <=> $b } keys %Map) {

	my $lgDir = '$dOut/LG$lg';	
	mkdir $lgDir unless (-d $lgDir);

	foreach my $sex (keys %{$Map{$lg}}) {
		
		my @marker = sort { $Map{$lg}{$sex}{'order'}{$a} <=> $Map{$lg}{$sex}{'order'}{$b} } keys %{$Map{$lg}{$sex}{'order'}};

		## basic statistics 

		$EvaluationStat{$lg}{$sex}{'markerNum'} = @marker;
		$EvaluationStat{$lg}{$sex}{'totalCM'} = abs($Map{$lg}{$sex}{'cm'}{$marker[-1]} - $Map{$lg}{$sex}{'cm'}{$marker[0]});
		$EvaluationStat{$lg}{$sex}{'averCM'} = sprintf("%.2f",$EvaluationStat{$lg}{$sex}{'totalCM'} / (@marker - 1));
		$EvaluationStat{$lg}{$sex}{'gapltx'} = gap_ltx_ratio(\%{$Map{$lg}{$sex}{'cm'}},$gapThreshold);

		## evaluate entire map quality
		@{$EvaluationStat{$lg}{$sex}{'chisqure'}} = calculate_chisquare(\%{$pwd{$lg}{$sex}},\%{$Map{$lg}{$sex}}); 
		$EvaluationStat{$lg}{$sex}{'RMSD'} = calculate_RMSD(\%{$Map{$lg}{$sex}},\%{$pwd{$lg}{$sex}}); 
		$EvaluationStat{$lg}{$sex}{'AFE'} = calculate_AFE(\%{$Map{$lg}{$sex}},\%{$pwd{$lg}{$sex}});
		
		print OUT ">",join("\t",("LG$lg",$sex,$EvaluationStat{$lg}{$sex}{'markerNum'},sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'totalCM'}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'averCM'}),sprintf("%.2f",$EvaluationStat{$lg}{$sex}{'gapltx'}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'chisqure'}->[0]),$EvaluationStat{$lg}{$sex}{'chisqure'}->[1],sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'chisqure'}->[2]),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'RMSD'}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'AFE'}))),"\n";

		## evaluate error of individual loci

		foreach my $marker (@marker) {

#			print $marker,"\n";

			## N.N.Fit (assess wheather the individual loci is properly placed)

			$EvaluationStat{$lg}{$sex}{'NNFit'}{$marker} = calculate_N_N_Fit(\%{$Map{$lg}{$sex}},\%{$pwd{$lg}{$sex}},$marker);

			## r_synteny --- its function is similar to N.N.Fit

			$EvaluationStat{$lg}{$sex}{'r_synteny'}{$marker} = calculate_r_synteny(\%{$Map{$lg}{$sex}},\%{$pwd{$lg}{$sex}},$marker);

			## mean chiSquare contribution --- large mcc symbolise the loci's influence to the entire map given the order of rest ones (cancel)

			## jump (standard for removing loci)--- large jump strongly accommend a removal of current loci

			$EvaluationStat{$lg}{$sex}{'jump'}{$marker} = calculate_jump(\%{$Map{$lg}{$sex}},\%{$pwd{$lg}{$sex}},$marker,\@{$EvaluationStat{$lg}{$sex}{'chisqure'}});

			## error Geno.prob --- when above index stay at a resonable interval,the false negative and false positive rate is decreased simutanouesly when perform imputation of missing data  or correction of error genotypes 
			
			print OUT join("\t",("LG$lg",$sex,$Map{$lg}{$sex}{'order'}{$marker},sprintf("%20s",$marker),$loc{$lg}{$sex}{$marker}{'type'},sprintf("%.3f",$Map{$lg}{$sex}{'cm'}{$marker}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'NNFit'}{$marker}{'r_unit'}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'NNFit'}{$marker}{'cm_unit'}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'r_synteny'}{$marker}),sprintf("%.3f",$EvaluationStat{$lg}{$sex}{'jump'}{$marker}))),"\n";
		
		}

	}
	
}
close (OUT) ;
store \%EvaluationStat,"$dOut/$fKey.eMap.hash"; # store the data structrure of statistic in a file 
#------------------------------------------------------------------
# draw 
#------------------------------------------------------------------

print "draw ...\n";
## entire quality_ RMSD

mkdir "$dOut/RMSD" unless (-d "$dOut/RMSD") ;
open (OUT,">$dOut/RMSD/$fKey.RMSD.list") or die $!;
{	## 全局参数
	my ($xMax) = sort {$b <=> $a} keys %Map;

	my ($yMax) = sort {$b <=> $a} ((map {$EvaluationStat{$_}{'sexAver'}{'RMSD'}} keys %Map),(map {$EvaluationStat{$_}{'female'}{'RMSD'}} keys %Map),(map {$EvaluationStat{$_}{'male'}{'RMSD'}} keys %Map));
	$yMax += 5;
	my $yStep = max(int(($yMax/5) - ($yMax/5)%5),5);
		
	print OUT "Type:Simple\n";                                                   
	print OUT "Width:800\n";
	print OUT "Height:600\n";
	print OUT "FontSize:32\n";
	print OUT "WholeScale:0.9\n";
	print OUT "ScaleLen:8\n";
	print OUT "OffsetPer:0.25\n";
	print OUT "UnitPer:0.25\n";
	print OUT "MovePer:0.25\n";
	print OUT "XScalePos:0.5\n";
	print OUT "XScaleRoate:75\n";
	print OUT "XStart:0\n";
	print OUT "XStep:1\n";
	print OUT "XEnd:$xMax\n";
	print OUT "YStart:0\n";
	print OUT "YEnd:$yMax\n";
	print OUT "YStep:$yStep\n";
	print OUT "XUnit:1\n";
	print OUT "MarkPos:rt\n";
	print OUT "MarkScale:0.6\n";
	print OUT "Note:Root mean square difference\n";
	print OUT "X:\n";
	print OUT "Y: # RMSD\n";
	print OUT "Scale:\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		print OUT "LG$lg\n";
	}
	print OUT ":End\n";
	print OUT "\n\n";

	print OUT "Color:#2292DD\n";
	print OUT "Mark:female\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:$EvaluationStat{$lg}{'female'}{'RMSD'}\n";
	}

	print OUT "\n\n";

	print OUT "Color:#CC3333\n";
	print OUT "Mark:sexAver\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:$EvaluationStat{$lg}{'sexAver'}{'RMSD'}\n";
	}
	print OUT "\n\n";

	print OUT "Color:#52CC33\n";
	print OUT "Mark:male\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:$EvaluationStat{$lg}{'male'}{'RMSD'}\n";
	}

	print OUT "\n\n";
	close (OUT);

	my $cur_dir = `pwd`;chomp $cur_dir;
	my $dirname = Cwd::abs_path("$dOut/RMSD");
		
	chdir $dirname;
	`perl $Bin/distributing_svg.pl "$fKey.RMSD.list" "$fKey.RMSD.list.svg"`;
	`perl $Bin/svg2xxx_release/svg2xxx "$fKey.RMSD.list.svg"`;

	chdir $cur_dir;
}
close (OUT) ;

## entire quality_ AFE

mkdir "$dOut/AFE" unless (-d "$dOut/AFE") ;
open (OUT,">$dOut/AFE/$fKey.AFE.list") or die $!;
{	## 全局参数
	my ($xMax) = sort {$b <=> $a} keys %Map;
		
	print OUT "Type:Simple\n";                                                   
	print OUT "Width:800\n";
	print OUT "Height:600\n";
	print OUT "FontSize:32\n";
	print OUT "WholeScale:0.9\n";
#	print OUT "MultiY:1\n";
	print OUT "ScaleLen:8\n";
	print OUT "OffsetPer:0.25\n";
	print OUT "UnitPer:0.25\n";
	print OUT "MovePer:0.25\n";
	print OUT "XScalePos:0.5\n";
	print OUT "XScaleRoate:75\n";
	print OUT "XStart:0\n";
	print OUT "XStep:1\n";
	print OUT "XEnd:$xMax\n";
	print OUT "YStart:0\n";
	print OUT "YEnd:100\n";
	print OUT "YStep:20\n";
	print OUT "XUnit:1\n";
	print OUT "MarkPos:rt\n";
	print OUT "MarkScale:0.6\n";
	print OUT "Note:Average fractional error\n";
	print OUT "X:\n";
	print OUT "Y: % AFE\n";
	print OUT "Scale:\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		print OUT "LG$lg\n";
	}
	print OUT ":End\n";
	print OUT "\n\n";

	print OUT "Color:#2292DD\n";
	print OUT "Mark:female\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:",100*$EvaluationStat{$lg}{'female'}{'AFE'},"\n";
	}

	print OUT "\n\n";

	print OUT "Color:#CC3333\n";
	print OUT "Mark:sexAver\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:",100*$EvaluationStat{$lg}{'sexAver'}{'AFE'},"\n";
	}
	print OUT "\n\n";

	print OUT "Color:#52CC33\n";
	print OUT "Mark:male\n";

	foreach my $lg (sort {$a <=> $b} keys %Map) {
		my $x = $lg - 1;
		print OUT "$x:",100*$EvaluationStat{$lg}{'male'}{'AFE'},"\n";
	}

	print OUT "\n\n";
	close (OUT);

	my $cur_dir = `pwd`;chomp $cur_dir;
	my $dirname = Cwd::abs_path("$dOut/AFE");
		
	chdir $dirname;
	`perl $Bin/distributing_svg.pl "$fKey.AFE.list" "$fKey.AFE.list.svg"`;
	`perl $Bin/svg2xxx_release/svg2xxx "$fKey.AFE.list.svg"`;

	chdir $cur_dir;
}
close (OUT) ;

## poor_fit marker detection

foreach my $lg (sort { $a <=> $b } keys %EvaluationStat) {

	my $lgDir = "$dOut/LG$lg";
	mkdir $lgDir unless (-d $lgDir);	
	mkdir "$lgDir/Result" unless (-d "$lgDir/Result");
	mkdir "$lgDir/Result/poorFit" unless (-d "$lgDir/Result/poorFit");


	foreach my $sex (keys %{$EvaluationStat{$lg}}) {

		my $poor_fit = locatePoorFit(\%{$Map{$lg}{$sex}},\%{$EvaluationStat{$lg}{$sex}},$sex); 

		if (scalar keys %{$poor_fit} > 0) {

			open (OUT,">$lgDir/Result/poorFit/LG$lg.$sex.poorFit.info") or die $!;
			foreach my $marker (sort {$poor_fit->{$b} <=> $poor_fit->{$a}} keys %{$poor_fit}) {
				print OUT "#\tLG$lg\t",$sex,"\t",$marker,"\t",sprintf("%.3f",$poor_fit->{$marker}),"\n";
			}

			my $mean = sum(values %{$poor_fit}) / (scalar keys %{$poor_fit});
			my @poor_fit = grep {$poor_fit->{$_} > $mean} keys %{$poor_fit};
			print OUT join("\t",("LG$lg",$sex,$_)),"\n" for(sort {$poor_fit->{$b} <=> $poor_fit->{$a}} @poor_fit);

			## seek ties (to be continued...)
			
			close (OUT);

		}

		## draw N.N.Fit histogram

		open (OUT,">$lgDir/Result/poorFit/LG$lg.$sex.poorFit.list") or die $!;

		## 全局参数
		my @marker = sort {$Map{$lg}{$sex}{'order'}{$a} <=> $Map{$lg}{$sex}{'order'}{$b}} keys %{$Map{$lg}{$sex}{'order'}};
		my $xMax = @marker;
		my ($yMax) = sort { $b <=> $a } map {$EvaluationStat{$lg}{$sex}{'NNFit'}{$_}{'cm_unit'}} @marker;
		$yMax += 5;
		
		my $yStep = max(int (($yMax/5)-($yMax/5)%5),5);

		my $width = 40 * @marker;
		my $height = 3*$width/4;

		print OUT "Type:Simple\n";                                                   
#		print OUT "PointSize:3\n";
		print OUT "Width:$width\n";
		print OUT "Height:$height\n";
		print OUT "FontSize:32\n";
		print OUT "WholeScale:0.9\n";
		print OUT "ScaleLen:8\n";
		print OUT "UnitPer:0.8\n";
		print OUT "MovePer:0.5\n";
		print OUT "XScalePos:0.5\n";
		print OUT "XScaleRoate:75\n";
		print OUT "XStart:0\n";
		print OUT "XStep:1\n";
		print OUT "XEnd:$xMax\n";
		print OUT "YStart:0\n";
		print OUT "YEnd:$yMax\n";
		print OUT "YStep:$yStep\n";
		print OUT "XUnit:1\n";
		print OUT "MarkPos:rt\n";
		print OUT "MarkScale:0.6\n";
		print OUT "Note:Distribution of Nearest neighbour fit\n";
		print OUT "X:\n";
		print OUT "Y: N.N.Fit(CM)\n";
		print OUT "Scale:\n";

		foreach my $marker (@marker) {

			print OUT "$marker\n";

		}
		print OUT ":End\n";
		print OUT "\n\n";

		print OUT "Color:#9B97C8\n";
		print OUT "Mark:LG$lg $sex\n";

		for (my $i=0;$i<@marker ;$i++) {

			print OUT "$i:$EvaluationStat{$lg}{$sex}{'NNFit'}{$marker[$i]}{'cm_unit'}\n";
		}

		print OUT "\n\n";
		close (OUT);

		my $cur_dir = `pwd`;chomp $cur_dir;
		my $dirname = Cwd::abs_path("$lgDir/Result/poorFit");
		
		chdir $dirname;
		`perl $Bin/distributing_svg.pl "LG$lg.$sex.poorFit.list" "LG$lg.$sex.poorFit.list.svg"`;
		`perl $Bin/svg2xxx_release/svg2xxx "LG$lg.$sex.poorFit.list.svg"`;

		chdir $cur_dir;
	}
	
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub locatePoorFit {#
	my ($ref_map,$ref_eMap,$sex) = @_;

	my @all_loci = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} keys %{$ref_map->{'order'}};
	my %result = ();

	foreach my $loci (@all_loci) {

		if ($ref_map->{'order'}{$loci} == 1) {
			
			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[1]}{'cm_unit'} > 0) {
				
				my $ratio = $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} / $ref_eMap->{'totalCM'};
				$result{$loci} = $ratio;
			}

		}elsif($ref_map->{'order'}{$loci} == @all_loci){

			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[@all_loci-1]}{'cm_unit'} > 0) {
				
				my $ratio = $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} / $ref_eMap->{'totalCM'};
				$result{$loci} = $ratio;
			}

		}else{

			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}-2]}{'cm_unit'} > 0
				&& $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}]}{'cm_unit'} > 0) {

				my $ratio = $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} / $ref_eMap->{'totalCM'};
				$result{$loci} = $ratio;
			}
		}
	}

	return \%result;
}

sub locatePoorFit_old {#
	my ($ref_map,$ref_eMap,$sex,$cm_threshold,$cv_threshold) = @_;

	my @all_loci = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} keys %{$ref_map->{'order'}};
	my %suspicious = map {$_,1} grep {$ref_eMap->{'NNFit'}{$_}{'cm_unit'} > $cm_threshold} @all_loci;

	my @result = ();

	foreach my $loci (keys %suspicious) {

		if ($ref_map->{'order'}{$loci} == 1) {
			
			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[1]}{'cm_unit'} > 0
				&& ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[1]}{'cm_unit'})/$ref_eMap->{'NNFit'}{$loci}{'cm_unit'} > $cv_threshold) {

				push @result,$loci;
			}

		}elsif($ref_map->{'order'}{$loci} == @all_loci){

			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[@all_loci-1]}{'cm_unit'} > 0
				&& ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[@all_loci-1]}{'cm_unit'})/$ref_eMap->{'NNFit'}{$loci}{'cm_unit'} > $cv_threshold) {
				
				push @result,$loci;
			}

		}else{

			if ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}-2]}{'cm_unit'} > 0
				&& $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}]}{'cm_unit'} > 0
				&& ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}-2]}{'cm_unit'}) / $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} > $cv_threshold
				&& ($ref_eMap->{'NNFit'}{$loci}{'cm_unit'} - $ref_eMap->{'NNFit'}{$all_loci[$ref_map->{'order'}{$loci}]}{'cm_unit'}) / $ref_eMap->{'NNFit'}{$loci}{'cm_unit'} > $cv_threshold) {

				push @result,$loci;
			}
		}
	}

	return @result;
}
sub gap_ltx_ratio {
	my $ref_map = shift;
	my $theshold = shift;

	my @cm = sort { $a <=> $b } values %{$ref_map};
	my @gap_ltx = grep {$_ < $theshold} map {$cm[$_] - $cm[$_-1]} 1..@cm-1;

	my $result = sprintf("%.4f",@gap_ltx / (@cm -1));
	return $result;

}

sub calculate_r_synteny {#
	my ($ref_map,$ref_pwd,$loci) = @_;
	
	my @valid = grep {if(defined $ref_pwd->{$loci}{$_}){$ref_pwd->{$loci}{$_}{'r'} < $rThreshold && $ref_pwd->{$loci}{$_}{'lod'} > $lodThreshold}} keys %{$ref_map->{'order'}};
	
	if (@valid != 0) {

		my @upStream = grep {$ref_map->{'order'}{$_} < $ref_map->{'order'}{$loci}} @valid;
		my @downStream = grep {$ref_map->{'order'}{$_} > $ref_map->{'order'}{$loci}} @valid;


		my $syntenyNum = 0;
		if (@upStream != 0) {


			my @mapOrder = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} @upStream;
			my @observedOrder = sort {$ref_pwd->{$loci}{$b}{'r'} <=> $ref_pwd->{$loci}{$a}{'r'}} @upStream;

			my @synteny = pickSyntenyLoci(\@mapOrder,\@observedOrder);

			$syntenyNum += @synteny;
		}

		if (@downStream != 0) {

			my @mapOrder = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} @downStream;
			my @observedOrder = sort {$ref_pwd->{$loci}{$a}{'r'} <=> $ref_pwd->{$loci}{$b}{'r'}} @downStream;

			my @synteny = pickSyntenyLoci(\@mapOrder,\@observedOrder);

			$syntenyNum += @synteny;

		}

		return $syntenyNum / @valid;


	}else{
		return 0;
	}
}

sub calculate_jump {#
	
	my ($ref_map,$ref_pwd,$loci,$chi_square) = @_;
	
	my $cm = $ref_map->{'cm'}{$loci};
	my $order = $ref_map->{'order'}{$loci};
	
	delete $ref_map->{'cm'}{$loci};
	delete $ref_map->{'order'}{$loci};

	my $jump = 0;
	my @chi_square = calculate_chisquare($ref_pwd,\%{$ref_map});

	$ref_map->{'cm'}{$loci} = $cm;
	$ref_map->{'order'}{$loci} = $order;
	
	if ($chi_square->[1] < $chi_square[1]) {
		print "insufficient linkage in data\n";
		exit(-1);
	}

	$jump = ($chi_square->[0] - $chi_square[0] + $chi_square[1] - $chi_square->[1])/sqrt(2*($chi_square->[1] - $chi_square[1]));
	return $jump ;
}

sub PickMaxGroup {
    my(@Pos)=@_;
    
    my @dp;#[$pos][0]-total [$pos][1]-prePos
    my $MaxEnd=0;
    my $MaxScore=0;
    for(my $i=0;$i<@Pos;$i++)
    {
        my $max_total=1;
        my $pre=-1;
        for(my $j=$i-1;$j>=0;$j--)
        {
            if($Pos[$i]>=$Pos[$j])
            {
                if($dp[$j][0]>=$max_total)
                {
                    $max_total=1+$dp[$j][0];
                    $pre=$j;
                    #last;
                }
            }
        }
        $dp[$i][0]=$max_total;
        $dp[$i][1]=$pre;
        if($MaxScore<$max_total)
        {
            $MaxScore=$max_total;
            $MaxEnd=$i;
        }
    }
    
    #find max
    my @result;
    for(my $i=$MaxEnd;$i>=0;)
    {
        push @result,$i;
        $i=$dp[$i][1];
    }
    return reverse @result;
}

sub pickSyntenyLoci {#
	my ($baseline,$arr)=@_;
	my @originalOrder=@{$arr};
	my @reversedOrder=reverse @{$arr};

	## original order

	my @posInBase=map {my $x=$_;grep {$baseline->[$_] eq $originalOrder[$x]}0..$#originalOrder;}0..$#originalOrder;

	my @originalSynteny=PickMaxGroup(@posInBase);

	## reverse order
	
	my @reversePosInBase=map {my $x=$_;grep {$baseline->[$_] eq $reversedOrder[$x]}0..$#reversedOrder;}0..$#reversedOrder;

	my @reverseSynteny=PickMaxGroup(@reversePosInBase);

	##compare

	my @result=();

	if (@originalSynteny >= @reverseSynteny) {
		@result = map {$baseline->[$posInBase[$_]]} @originalSynteny;
	}else{
		@result = map {$baseline->[$reversePosInBase[$_]]} @reverseSynteny;
	}
	return @result;
}

sub calculate_N_N_Fit {#
	my ($ref_map,$ref_pwd,$loci) = @_;

	my %result = ();

	if ($ref_map->{'order'}{$loci} == 1) {

		foreach my $marker (sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} keys %{$ref_map->{'order'}}) {

			if (exists $ref_pwd->{$loci}{$marker}) {

				$result{'r_unit'} = abs($ref_pwd->{$loci}{$marker}{'r'} - &{$inverseMapFunction{$mapFunction}}($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
				$result{'cm_unit'} = abs(&{$mapFunction{$mapFunction}}($ref_pwd->{$loci}{$marker}{'r'}) - abs($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));

				last;
			}
		}
		
	}elsif($ref_map->{'order'}{$loci} == scalar keys %{$ref_map->{'order'}}){

		foreach my $marker (sort {$ref_map->{'order'}{$b} <=> $ref_map->{'order'}{$a}} keys %{$ref_map->{'order'}}) {

			if (exists $ref_pwd->{$loci}{$marker}) {

				$result{'r_unit'} = abs($ref_pwd->{$loci}{$marker}{'r'} - &{$inverseMapFunction{$mapFunction}}($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
				$result{'cm_unit'} = abs(&{$mapFunction{$mapFunction}}($ref_pwd->{$loci}{$marker}{'r'}) - abs($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
				
				last;
			}
		}
	}else{

		my @downStream = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} grep {$ref_map->{'order'}{$_} > $ref_map->{'order'}{$loci}} keys %{$ref_map->{'order'}};
		my @upStream = sort {$ref_map->{'order'}{$b} <=> $ref_map->{'order'}{$a}} grep {$ref_map->{'order'}{$_} < $ref_map->{'order'}{$loci}} keys %{$ref_map->{'order'}};
		
		foreach my $marker (@downStream) {

			next unless (exists $ref_pwd->{$loci}{$marker}) ;
			
			$result{'r_unit'} += abs($ref_pwd->{$loci}{$marker}{'r'} - &{$inverseMapFunction{$mapFunction}}($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
			$result{'cm_unit'} += abs(&{$mapFunction{$mapFunction}}($ref_pwd->{$loci}{$marker}{'r'}) - abs($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
			last;

		}

		foreach my $marker (@upStream) {

			next unless (exists $ref_pwd->{$loci}{$marker}) ;
			
			$result{'r_unit'} += abs($ref_pwd->{$loci}{$marker}{'r'} - &{$inverseMapFunction{$mapFunction}}($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
			$result{'cm_unit'} += abs(&{$mapFunction{$mapFunction}}($ref_pwd->{$loci}{$marker}{'r'}) - abs($ref_map->{'cm'}{$marker} - $ref_map->{'cm'}{$loci}));
			
			last;
		}
	
	}

	return \%result;

}

sub calculate_AFE {
	my ($ref_map,$ref_pwd) = @_;

	my $total_obs_dist = 0;
	my $total_diff_dist = 0;

	my @marker = sort { $ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b} } keys %{$ref_map->{'order'}};

	for (my $i=0;$i<@marker-1 ;$i++) {
		for (my $j=$i+1;$j<@marker ;$j++) {

			next unless (defined $ref_pwd->{$marker[$i]}{$marker[$j]});

			$total_obs_dist += &{$mapFunction{$mapFunction}}($ref_pwd->{$marker[$i]}{$marker[$j]}{'r'});
			$total_diff_dist += abs(&{$mapFunction{$mapFunction}}($ref_pwd->{$marker[$i]}{$marker[$j]}{'r'}) - $ref_map->{'cm'}{$marker[$j]} + $ref_map->{'cm'}{$marker[$i]});
			
		}
	}

	my $afe = $total_diff_dist / $total_obs_dist;
	return $afe;
}
sub calculate_RMSD {
	my ($ref_map,$ref_pwd) = @_;

	my $pair_num = 0;
	my $rmsd = 0;

	my @marker = sort {$ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b}} keys %{$ref_map->{'order'}};

	for (my $i=0;$i<@marker-1 ;$i++) {
		for (my $j=$i+1;$j<@marker ;$j++) {
#			next unless (defined $ref_pwd->{$marker[$i]}{$marker[$j]});
			if (defined $ref_pwd->{$marker[$i]}{$marker[$j]}) {

				$pair_num++;
				$rmsd += (&{$mapFunction{$mapFunction}}($ref_pwd->{$marker[$i]}{$marker[$j]}{'r'}) - $ref_map->{'cm'}{$marker[$j]} + $ref_map->{'cm'}{$marker[$i]})**2;
			}
			
			
		}
	}

	$rmsd = sqrt($rmsd / $pair_num);

	return $rmsd;

}
sub calculate_chisquare{# calculate the value of  goodness-of-fit statistic
	
	my ($pairWise,$ref_map)=@_;
	my $G_square=0;

	my $order = [sort { $ref_map->{'order'}{$a} <=> $ref_map->{'order'}{$b} } keys %{$ref_map->{'order'}}];
	my $distance = [map {$ref_map->{'cm'}{$order->[$_]} - $ref_map->{'cm'}{$order->[$_-1]}} 1..@{$order}-1];
	
	if (scalar @{$order}<3) {
		warn "Error,the marker list must contain more than three markers\n";
		exit(-1);
	}
	
	my $df=&calculateValidPairs(\%{$pairWise},\@{$order})-@{$order}+1;
	return (0,0,0) if ($df == 0) ;
	
	for (my $i=0;$i<@{$order}-1 ;$i++) {
		for (my $j=$i+1;$j<@{$order} ;$j++) {
			
			next if (!defined $pairWise->{$order->[$i]}{$order->[$j]}) ;
			
			if ($pairWise->{$order->[$i]}{$order->[$j]}{"r"} <= $rThreshold and $pairWise->{$order->[$i]}{$order->[$j]}{"lod"} >= $lodThreshold) {
				
				if ($pairWise->{$order->[$i]}{$order->[$j]}{"r"} == 0) {
					
					my $Nt=$pairWise->{$order->[$i]}{$order->[$j]}{"lod"}/log10(2);
					my $Nr=0;
					my $Nn=$Nt;

					my $mapDist = sum(@{$distance}[$i..$j-1]);

					$mapDist = 0.000001 if ($mapDist == 0) ;

					$G_square+=2*(0-$Nn*log(1-&{$inverseMapFunction{$mapFunction}}($mapDist)));
				
				}else{
					
					my $Nt=$pairWise->{$order->[$i]}{$order->[$j]}{"lod"}/(log10(2)+$pairWise->{$order->[$i]}{$order->[$j]}{"r"}*log10($pairWise->{$order->[$i]}{$order->[$j]}{"r"})+(1-$pairWise->{$order->[$i]}{$order->[$j]}{"r"})*log10(1-$pairWise->{$order->[$i]}{$order->[$j]}{"r"}));
					my $Nr=$Nt*$pairWise->{$order->[$i]}{$order->[$j]}{"r"} ;
					my $Nn=$Nt*(1-$pairWise->{$order->[$i]}{$order->[$j]}{"r"});
					#The above three scalars  represent virtual number of total gametes,recombinatorial gametes,non_recombinatorial gametes respectively and sequentialy calculated by recombination frequency and LOD score between pairs of loci
					

					my $mapDist = sum(@{$distance}[$i..$j-1]);
					$mapDist = 0.000001 if ($mapDist == 0) ;

					if ($mapDist < 0) {

						print $mapDist,"\t",$i,"\t",$j,"\t",$order->[$i],"\t",$order->[$j],"\n";
						print "error,negative distance during caculate chisquare\n";
						exit(-1);
					}
#					print join("\t",($pairWise->{$order->[$i]}{$order->[$j]}{"r"},&{$inverseMapFunction{$mapFunction}}(sum(@{$distance}[$i..$j-1])),1-&{$inverseMapFunction{$mapFunction}}(sum(@{$distance}[$i..$j-1])))),"\n";
					$G_square+=2*($Nr*log($pairWise->{$order->[$i]}{$order->[$j]}{"r"})+$Nn*log(1-$pairWise->{$order->[$i]}{$order->[$j]}{"r"})-$Nr*log(&{$inverseMapFunction{$mapFunction}}($mapDist))-$Nn*log(1-&{$inverseMapFunction{$mapFunction}}($mapDist))) ;
				}	
			}
		}
	}

	my $mean = $G_square / $df;
	
	return ($G_square,$df,$mean) ;
}

sub log10 {#
	my $l=shift;
	return log($l)/log(10);
}

sub sum {#
	my $sum=0;
	$sum+=$_  foreach (@_) ;
	return $sum;
}

sub calculateValidPairs {#
	
	my ($pairWise,$markerList)=@_;
	my $valid=0;
	
	for (my $i=0;$i<@{$markerList}-1 ;$i++) {
		for (my $j=$i+1;$j<@{$markerList} ;$j++) {
			
			next if (!defined $pairWise->{$markerList->[$i]}{$markerList->[$j]}) ;
			$valid++ if ($pairWise->{$markerList->[$i]}{$markerList->[$j]}{"r"} <= $rThreshold and $pairWise->{$markerList->[$i]}{$markerList->[$j]}{"lod"} >= $lodThreshold) ;
		
		}
	}
	
	return $valid;
}

sub loadMap {#
	my ($refMap,$fIn)=@_;
	if (-d $fIn) {
		my @mapFile = glob("$fIn/*.map");
		foreach my $map (@mapFile) {
			my ($lgID,$sex)=basename($map)=~/^\D+(\d+)\.(\w+)\.map$/;
			open (IN,"$map") or die $!;
			$/="\n";
			my $order=0;
			while (<IN>) {
				chomp;
				s/\r//g;
				next if (/^$/ || /^;/ || /^group/) ;
				my ($marker,$cm) = split;
				$refMap->{$lgID}{$sex}{'cm'}{$marker}=$cm;
				$refMap->{$lgID}{$sex}{'order'}{$marker}=++$order;
			}
			close (IN) ;
		}
	}else{
		print "your input should be a directory,please check!\n";
		exit(-1);
	}

	# check wheather the map is complete
	foreach my $lg (keys %{$refMap}) {
		next if (scalar keys %{$refMap->{$lg}} == 3) ;
		print "LG$lg map:",join("\t",keys %{$refMap->{$lg}}),"\n";
		exit(-1);
	}
}

sub loadPWD {#
	my ($refPwd,$fIn)=@_;
	if (-d $fIn) {
		my @file = glob("$fIn/*.pwd");
		foreach my $pwd (@file) {
			my ($groupNo,$sex)=basename($pwd)=~/^\D+(\d+)\.(\w+)\.pwd$/;
			open (IN,"$pwd") or die $!;
			$/="\n";
			while (<IN>) {
				chomp;
				next if (/^$/ or /^;/) ;
				next unless (/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) ;
				$refPwd->{$groupNo}{$sex}{$1}{$2}{"r"}=$3;
				$refPwd->{$groupNo}{$sex}{$2}{$1}{"r"}=$3;
				$refPwd->{$groupNo}{$sex}{$1}{$2}{"lod"}=$4;
				$refPwd->{$groupNo}{$sex}{$2}{$1}{"lod"}=$4;
			}
			close (IN) ;
		}
	}else{
		print "your input should be a directory,please check!\n";
		exit(-1);
	}

	 
	# check integrity of pwd file

	foreach my $lg (keys %{$refPwd}) {

		next if (scalar keys %{$refPwd->{$lg}} == 3) ;
		print "LG$lg pwd:",join("\t",keys %{$refPwd->{$lg}}),"\n";
		exit(-1);

	}
	
}


sub loadLOC {#
	my ($refLoc,$fIn)=@_;
	if (-d $fIn) {
		my @file = glob("$fIn/*.loc");
		foreach my $loc (@file) {
			my ($lg,$sex) = basename($loc) =~/\D+(\d+)\.(\w+)\.loc$/;
			open (IN,"$loc") or die $!;
			$/="\n";
			while (<IN>) {
				chomp;
				next if (/^;/ or /^$/ or /^name|^popt|^nind|^nloc/) ;
				my ($marker,$crosstype,$lp,@geno)=split;
				$crosstype=~s/<|>//g;
				$lp=~s/{|}//g;
				$refLoc->{$lg}{$sex}{$marker}{'type'}=$crosstype;
				$refLoc->{$lg}{$sex}{$marker}{'lp'}=$lp;
				$refLoc->{$lg}{$sex}{$marker}{'proGeno'}=\@geno;
			}
			close (IN) ;
		}
	}else{
		print "your input should be a directory,please check!\n";
		exit(-1);
	}

	 
	# check integrity of loc file

	foreach my $lg (keys %{$refLoc}) {

		next if (scalar keys %{$refLoc->{$lg}} == 3) ;
		print "LG$lg loc:",join("\t",keys %{$refLoc->{$lg}}),"\n";
		exit(-1);

	}
}

sub Haldane {#
	my $r=shift;
	my $result=-100*(1/2)*log(1-2*$r);
	return $result;
}

sub Kosambi {#
	my $r=shift;
	my $result=100*(1/4)*log((1+2*$r)/(1-2*$r));
	return $result;
}

sub inverseHaldane {#
	my $r=shift;
	my $result=(1/2)*(1-exp((-1)*$r/50));
	return $result;
}

sub inverseKosambi {#
	my $r=shift;
	my $result=(1/2)*(exp($r/25)-1)/(exp($r/25)+1);
	return $result;
}

sub max{
	my ($a,$b) = @_;
	return ($a<$b)?$b:$a;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $0 
Version: $version
Contact: Ma Chouxian <macx\@biomarker.com.cn> 
Description:
	This program,recieves a directory of linakge data as input,is used to evaluate the liankage genetic maps of CP population.
	
	Files in the input directory include map,loc and pwd,named like :LG+.{sexAver,male,female}.{pwd,loc,map}.

	A thourough statistic of map is performed by using several basic statistics such as the number of loci,total and average map distance of linkage group,
	ratio of gap less than given centiMorgan threshold,and some sophisticated statistics ie. final chisquare,RMSD,AFE are used to evaluate the general map quality.
			 
	To evaluate the location of individual loci, a statistic named N.N.Fit is introduced from joinmap to detect poor fit loci in the final map of corresbonding linkage group.

	Distribution of RMSD,AFE,N.N.Fit are illustrated in a graphic way.(see output directory/RMSD,output directory/AFE and output directory/LG+/Result/poorFit/*.png respectively)

	A list of suspicious markers are provided finally(see output directory/LG+/Result/poorFit/*.info).It\'s a pity that this program fails to detect poor fit regions of input maps.

Usage:
  Options:
  -i		<dir>		Directory of linkage analysis data, forced
  -k		<str>		Key of output file,forced
  -d		<str>		Directory where output file produced,optional,default [./]
  
  -gap		<int>		Gap threshold of map,default [5]
  -mapFunction	<str>		Map function used,default ["Kosambi"]
  -r		<float>		Maximum r threshold,default [0.4]
  -lod		<float>		Minimum lod threshold,default [1]
  

  -h         Help

USAGE
	print $usage;
	exit;
}


