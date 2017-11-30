#!/usr/Bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$Bin/Statistics-RankCorrelation-0.1205/lib/Statistics/";
use RankCorrelation;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($indir, $key, $fgenotype, $outdir, $noise_ratio, $process, $queue, $lineMapFlag, $mapEvaFlag);
GetOptions(
				"help|?" =>\&USAGE,
				"id:s"=>\$indir,
				"igeno:s"=>\$fgenotype,
				"k:s"=>\$key,
				"od:s"=>\$outdir,
				"r:f"=>\$noise_ratio,
				"mapEva:s"=>\$mapEvaFlag,

				"p:i"=>\$process,
				"q:i"=>\$queue,
			) or &USAGE;
&USAGE unless ($indir and $key);
# ------------------------------------------------------------------
# global 
# ------------------------------------------------------------------
$outdir||="$key"."Result";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$indir=AbsolutePath("dir",$indir);
$fgenotype=AbsolutePath("file",$fgenotype);

$noise_ratio||=0.05;
$process||=50;
$queue||="general.q";

my $workdir=`pwd`;
chomp $workdir;
my %lg;
# ------------------------------------------------------------------
# deal with
# ------------------------------------------------------------------
########## draw
my @file=glob("$indir/*.sexAver.map");
foreach my $file(@file) {
	my ($lg) = basename($file) =~/(\S+)\.sexAver\.map/;
	$lg{$lg} = 1;
}


open (SH,">$outdir/$key.draw.sh") or die $!;
foreach my $lg (keys %lg) {
	`mkdir $outdir/$lg` unless(-d "$outdir/$lg");
	`mkdir $outdir/$lg/Data` unless (-d "$outdir/$lg/Data");
	my @file=glob("$indir/$lg.*");
	foreach my $file (@file) {
		`cp $file $outdir/$lg/Data/`;
	}
	`mkdir $outdir/$lg/Result` unless (-d "$outdir/$lg/Result");
	`mkdir $outdir/$lg/Result/synteny` unless(-d  "$outdir/$lg/Result/synteny");
	
	`mkdir $outdir/$lg/Result/haploSource` unless (-d "$outdir/$lg/Result/haploSource") ;
	`mkdir $outdir/$lg/Result/heatMap` unless(-d "$outdir/$lg/Result/heatMap");
	######sexAver 

	print SH "perl $Bin/IndividualLinkagePhaseNumber.pl -i $indir/$lg.sexAver.loc -m $indir/$lg.sexAver.map -o $outdir/$lg/Result/haploSource/$lg.sexAver && \t";
	print SH "perl $Bin/statHaploSourceNoise.pl -i $outdir/$lg/Result/haploSource/$lg.sexAver.phase -od $outdir/$lg/Result/haploSource -r $noise_ratio&& \n";

	print SH "perl $Bin/drawHeatmap.pl -i $indir/$lg.sexAver.map -p $indir/$lg.sexAver.pwd -k $lg.sexAver -d $outdir/$lg/Result/heatMap &&\n";

	if ( (-f "$indir/$lg.female.map") && (-f "$indir/$lg.female.pwd") && (-f "$indir/$lg.female.loc")) {
		print SH "perl $Bin/IndividualLinkagePhaseNumber.pl -i $indir/$lg.female.loc -m $indir/$lg.female.map -o $outdir/$lg/Result/haploSource/$lg.female && \t";
		print SH "perl $Bin/statHaploSourceNoise.pl -i $outdir/$lg/Result/haploSource/$lg.female.phase -od $outdir/$lg/Result/haploSource -r $noise_ratio&& \n";	
		print SH "perl $Bin/drawHeatmap.pl -i $indir/$lg.female.map -p $indir/$lg.female.pwd -k $lg.female -d $outdir/$lg/Result/heatMap &&\n";
	}
	if (-f "$indir/$lg.male.map" && -f "$indir/$lg.male.pwd" && -f "$indir/$lg.male.loc") {
		print SH "perl $Bin/IndividualLinkagePhaseNumber.pl -i $indir/$lg.male.loc -m $indir/$lg.male.map -o $outdir/$lg/Result/haploSource/$lg.male && \t";
		print SH "perl $Bin/statHaploSourceNoise.pl -i $outdir/$lg/Result/haploSource/$lg.male.phase -od $outdir/$lg/Result/haploSource -r $noise_ratio&& \n";
		print SH "perl $Bin/drawHeatmap.pl -i $indir/$lg.male.map -p $indir/$lg.male.pwd -k $lg.male -d $outdir/$lg/Result/heatMap &&\n";
	}
	if ((-f "$indir/$lg.female.map") && (-f "$indir/$lg.male.map")) {
		print SH "perl $Bin/drawGT3Map.pl -i1 $indir/$lg.female.map -i2 $indir/$lg.sexAver.map -i3 $indir/$lg.male.map -o $lg -d $outdir/$lg/Result/synteny -g $fgenotype -s run && \n";
	}
	
}
close (SH) ;

print "sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh -maxproc $process -reqsub -queue $queue $outdir/$key.draw.sh\n";
print `sh /share/nas2/genome/bmksoft/tool/qsub_sge_plus/v1.0/qsub_sge.plus.sh -maxproc $process -reqsub -queue $queue $outdir/$key.draw.sh`;


# ------------------------------------------------------------------
# summarise
# ------------------------------------------------------------------
my %haplo;
my %relation;
my %synteny;
my %info;
my %Map;
my %spearman_coef;
foreach my $lg (keys %lg) {
	###haploSource
	my @fhaplo_stat = glob ("$outdir/$lg/Result/haploSource/$lg.*.phase.noiseStat");
	foreach my $fhaplo_stat (@fhaplo_stat) {
		my $fhaplo_stat_name = basename($fhaplo_stat);
		my ($group,$sex) = $fhaplo_stat_name =~/(\w+)\.(\w+)\.phase\.noiseStat/;
		open (IN,$fhaplo_stat) or die $!;
		$/="\n";
		while (<IN>) {
			chomp;
			if(/Suspicious marker sum/){
				my ($nouse,$sum)=split /:/,$_;
				$haplo{$group}{$sex}{"susMarkerSum"}=$sum;
			}elsif (/gene sum/) {
				my ($nouse,$sum)=split /:/,$_;
				$haplo{$group}{$sex}{"geneSum"}=$sum;
			}elsif($_=~/noise sum/){
				my ($nouse,$sum)=split /:/,$_;
				$haplo{$group}{$sex}{"noiseSum"}=$sum;
			}elsif(/marker sum/){
				my ($nouse,$sum)=split /:/,$_;
				$info{$group}{$sex}{"markerSum"}=$sum;
			}elsif(/ind sum/){
				my ($nouse,$sum)=split /:/,$_;
				$info{$group}{$sex}{"indiSum"}=$sum;
			}
		}
		close (IN) ;
	}

	if (defined $lineMapFlag) {
		###relation
		my @frela_stat=glob("$outdir/$lg/Result/lineMap/$lg.*.relationStat");
		foreach my $frela_stat (@frela_stat) {
			my $frela_stat_name=basename($frela_stat);
			my ($group,$sex)=$frela_stat_name=~/(\w+)\.(\w+)\.relationStat/;
			open (IN,$frela_stat) or die $!;
			$/=">";
			while (<IN>) {
				chomp;
				next if(/^$/);
				my @line=split /\n/,$_;
				my $marker_sum=0;
				for (my $i=1;$i<@line;$i++) {
					my @line_marker=split /\t/,$line[$i];
					$marker_sum+=@line_marker;
				}
				$relation{$group}{$sex}=$marker_sum;
			}
			close (IN) ;
		}
	}

	###synteny
	my @fsynteny_stat=glob("$outdir/$lg/Result/synteny/*.synteny.stat.xls");
	foreach my $fsynteny_stat (@fsynteny_stat) {
		open (IN,$fsynteny_stat) or die $!;
		my $group;
		my $female_ratio;
		my $male_ratio;
		$/="\n";
		while (<IN>) {
			chomp;
			next if(/^$/);
			if (/[Ll][Gg](\w+)/i) {
				$group = $1;
			}elsif(/\//){
				my ($female_str,$male_str)=split /\s+/,$_;
				my @unit=split /\//,$female_str;
				if ($unit[1] ==0) {
					$female_ratio="--";
				}else{
					$female_ratio=$unit[0]/$unit[1];
				}

				@unit=split /\//,$male_str;
				if ($unit[1]==0) {
					$male_ratio="--";
				}else{
					$male_ratio=$unit[0]/$unit[1];
				}
				if ($group eq '') {
					print "$fsynteny_stat has some wrong!\n";
					exit (-1);
				}else{
					$synteny{$group}{"female"}=$female_ratio;
					$synteny{$group}{"male"}=$male_ratio;
				}
				$group='';
				$female_ratio=0;
				$male_ratio=0;
			}
		}
		close (IN) ;
	}

	### get map distance info 
	get_map_info("$indir/$lg.sexAver.map", \%{$Map{'sexAver'}}, $lg);
	get_map_info("$indir/$lg.male.map", \%{$Map{'male'}}, $lg);
	get_map_info("$indir/$lg.female.map", \%{$Map{'female'}}, $lg);

	### calculate spearman coefficient
	my @FS = map {$Map{'sexAver'}{$lg}{'order'}{$_}} sort {$Map{'female'}{$lg}{'order'}{$a} <=> $Map{'female'}{$lg}{'order'}{$b}} keys %{$Map{'female'}{$lg}{'order'}};
	my $c = Statistics::RankCorrelation->new( \@FS, [1..@FS-1]);
	my $n=$c->spearman;
	$spearman_coef{$lg}{'FS'} = abs($n);
	
	my @MS = map {$Map{'sexAver'}{$lg}{'order'}{$_}} sort {$Map{'male'}{$lg}{'order'}{$a} <=> $Map{'male'}{$lg}{'order'}{$b}} keys %{$Map{'male'}{$lg}{'order'}};
	$c = Statistics::RankCorrelation->new(\@MS, [1..@MS-1]);
	$n=$c->spearman;
	$spearman_coef{$lg}{'MS'} = abs($n);
}


# ------------------------------------------------------------------
# output
# ------------------------------------------------------------------
open (OUT,">$outdir/$key.quality.summarise") or die $!;

print OUT join("\t",
	"Group",
	"nind",
	"nloc_S",
	"nloc_F",
	"nloc_M",
	"md_S",
	"md_F",
	"md_M",
	"spearman_F",
	"spearman_M",
	"synteny_F",
	"synteny_M",
	"noise_S",
	"noise_F",
	"noise_M",
),"\n";
#if (defined $lineMapFlag) {
#	print OUT "\tinfo\t\t\t\tSynteny\t\tHaploSource\t\t\t\t\t\t\t\t\t\t\t\tRelation\t\n";
#	print OUT "\tall\tsexAver\tfemale\tmale\tfemale\tmale\tsexAver\t\t\t\tfemale\t\t\t\tmale\t\t\t\tsexAver\t\tfemale\t\tmale\t\t\n";
#	print OUT "group\tindi_sum\tmarker_sum\t\t\t";
#	print OUT "ratio\t\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\t";
#	print OUT "suspi_marker_sum\tratio\t";
#	print OUT "suspi_marker_sum\tratio\t";
#	print OUT "suspi_marker_sum\tratio\t\n";
#}else{
#	print OUT "\tinfo\t\t\t\tSynteny\t\tHaploSource\t\t\t\t\t\t\t\t\t\t\t\t\n";
#	print OUT "\tall\tsexAver\tfemale\tmale\tfemale\tmale\tsexAver\t\t\t\tfemale\t\t\t\tmale\t\t\t\t\n";
#	print OUT "group\tindi_sum\tmarker_sum\t\t\t";
#	print OUT "ratio\t\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\t";
#	print OUT "suspi_marker_sum\tgene_sum\tnoise_sum\tnoise_ratio\n";
#}
my @group=keys %lg;
my @groupSorted;
SortByNum("+",\@group,\@groupSorted);

for (my $i=0;$i<@groupSorted;$i++) {
	my $group=$groupSorted[$i];
	my $indi_sum='--';
	my $sexAver_marker_sum='--';
	my $female_marker_sum='--';
	my $male_marker_sum='--';

	my $female_map_distance = '--';
	my $male_map_distance = '--';
	my $sexAver_map_distance = '--';

	my $spearman_female = '--';
	my $spearman_male = '--';

	my $synteny_female='--';
	my $synteny_male='--';

	my $haplo_sexAver_susMarkerSum='--';
	my $haplo_sexAver_geneSum='--';
	my $haplo_sexAver_noiseSum='--';
	my $haplo_sexAver_noiseRatio='--';

	my $haplo_female_susMarkerSum='--';
	my $haplo_female_geneSum='--';
	my $haplo_female_noiseSum='--';
	my $haplo_female_noiseRatio='--';

	my $haplo_male_susMarkerSum='--';
	my $haplo_male_geneSum='--';
	my $haplo_male_noiseSum='--';
	my $haplo_male_noiseRatio='--';

	my $relation_sexAver='--';
	my $relation_sexAver_ratio='--';
	my $relation_female='--';
	my $relation_female_ratio='--';
	my $relation_male='--';
	my $relation_male_ratio='--';

	if (exists $haplo{$group}{"female"} && exists $haplo{$group}{"male"}) {
		$synteny_female=$synteny{$group}{"female"};
		$synteny_male=$synteny{$group}{"male"};

		$spearman_female = $spearman_coef{$group}{'FS'};
		$spearman_male = $spearman_coef{$group}{'MS'};
		$female_map_distance = $Map{'female'}{$group}{'cm'}[-1] - $Map{'female'}{$group}{'cm'}[0];
		$male_map_distance = $Map{'male'}{$group}{'cm'}[-1] - $Map{'male'}{$group}{'cm'}[0];
		$sexAver_map_distance = $Map{'sexAver'}{$group}{'cm'}[-1] - $Map{'sexAver'}{$group}{'cm'}[0];
	}

	if (exists $haplo{$group}{"female"}) {
		$female_marker_sum=$info{$group}{"female"}{"markerSum"};

		$haplo_female_susMarkerSum=$haplo{$group}{"female"}{"susMarkerSum"};
		$haplo_female_geneSum=$haplo{$group}{"female"}{"geneSum"};
		$haplo_female_noiseSum=$haplo{$group}{"female"}{"noiseSum"};

		$haplo_female_noiseRatio=$haplo{$group}{"female"}{"noiseSum"}/$haplo{$group}{"female"}{"geneSum"};
		if (defined $lineMapFlag) {
			$relation_female=$relation{$group}{"female"};
			$relation_female_ratio=$relation{$group}{"female"}/$info{$group}{"female"}{"markerSum"};
		}
	}
	if (exists $haplo{$group}{"male"}) {
		$male_marker_sum=$info{$group}{"male"}{"markerSum"};
		$haplo_male_susMarkerSum=$haplo{$group}{"male"}{"susMarkerSum"};
		$haplo_male_geneSum=$haplo{$group}{"male"}{"geneSum"};
		$haplo_male_noiseSum=$haplo{$group}{"male"}{"noiseSum"};
		$haplo_male_noiseRatio=$haplo{$group}{"male"}{"noiseSum"}/$haplo{$group}{"male"}{"geneSum"};

		if (defined $lineMapFlag) {
			$relation_male=$relation{$group}{"male"};
			$relation_male_ratio=$relation{$group}{"male"}/$info{$group}{"male"}{"markerSum"};
		}
	}

	if (exists $haplo{$group}{"sexAver"}{"susMarkerSum"}) {
		$indi_sum=$info{$group}{"sexAver"}{"indiSum"};
		$sexAver_marker_sum=$info{$group}{"sexAver"}{"markerSum"};

		$haplo_sexAver_susMarkerSum=$haplo{$group}{"sexAver"}{"susMarkerSum"};
		$haplo_sexAver_geneSum=$haplo{$group}{"sexAver"}{"geneSum"};
		$haplo_sexAver_noiseSum=$haplo{$group}{"sexAver"}{"noiseSum"};
		$haplo_sexAver_noiseRatio=$haplo{$group}{"sexAver"}{"noiseSum"}/$haplo{$group}{"sexAver"}{"geneSum"};

		if (defined $lineMapFlag) {
			$relation_sexAver=$relation{$group}{"sexAver"};
			$relation_sexAver_ratio=$relation{$group}{"sexAver"}/$info{$group}{"sexAver"}{"markerSum"};
		}
	}
	
	print OUT $group,"\t",$indi_sum,"\t",$sexAver_marker_sum,"\t",$female_marker_sum,"\t",$male_marker_sum,"\t";
	
	##map dist
	print OUT join("\t",
			$sexAver_map_distance,
			$female_map_distance,
			$male_map_distance,
		),"\t";

	## spearman coefficient && synteny ratio
	print OUT join("\t",
		sprintf("%.4f", $spearman_female),
		sprintf("%.4f", $spearman_male),
		sprintf("%.4f", $synteny_female),
		sprintf("%.4f", $synteny_male),
	),"\t";

	## noise
	print OUT join("\t",
			sprintf("%.6f", $haplo_sexAver_noiseRatio),
			sprintf("%.6f", $haplo_female_noiseRatio),
			sprintf("%.6f", $haplo_male_noiseRatio),
		),"\t";

#	print OUT $haplo_sexAver_susMarkerSum,"\t",$haplo_sexAver_geneSum,"\t",$haplo_sexAver_noiseSum,"\t",$haplo_sexAver_noiseRatio,"\t";
#	print OUT $haplo_female_susMarkerSum,"\t",$haplo_female_geneSum,"\t",$haplo_female_noiseSum,"\t",$haplo_female_noiseRatio,"\t";
#	print OUT $haplo_male_susMarkerSum,"\t",$haplo_male_geneSum,"\t",$haplo_male_noiseSum,"\t",$haplo_male_noiseRatio,"\t";	
	
	if (defined $lineMapFlag) {
		print OUT $relation_sexAver,"\t",$relation_sexAver_ratio,"\t";
		print OUT $relation_female,"\t",$relation_female_ratio,"\t";
		print OUT $relation_male,"\t",$relation_male_ratio,"\n";
	}else{
		print OUT "\n";
	}
}
close (OUT) ;


#######  thourough mapEvaluation
if (defined $mapEvaFlag) {
	print "perl $Bin/mapEvaluation.pl -i $indir -k $key -d $outdir...\n";
	`perl $Bin/mapEvaluation.pl -i $indir -k $key -d $outdir `;
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub get_map_info {#
	my ($fMap, $ref_map, $group) = @_;

	my $order = 1;
	open (IN,"$fMap") or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^;/) ;
		my ($marker, $cm) = split;
		push @{$ref_map->{$group}{'cm'}}, $cm;
		$ref_map->{$group}{'order'}{$marker} = $order++;
	}
	close (IN) ;

}
sub SortByNum {	###### 对那种以数字结尾的list按数字的大小排序, 排序类型用"+"表示升序，"-"表示降序。
				###### 例如：LG1，LG6，LG3，LG5 按升序排序后得到：LG1，LG3，LG5，LG6
	my ($sort_type,$ref_list,$ref_sorted)=@_;
	my %hash;
	my $sum=0;
	foreach my $name (@{$ref_list}) {
		$sum++;
		my ($word,$num)=$name=~/(\S+\D+)(\d+)/;
		$hash{$word}{$sum}=$num;
	}
	if ($sort_type eq "+") {
		foreach my $word (sort {$a cmp $b} keys %hash) {
			foreach my $sum (sort {$hash{$word}{$a} <=> $hash{$word}{$b}} keys %{$hash{$word}}) {
				push @{$ref_sorted},$word.$hash{$word}{$sum};
			}
		}
	}elsif($sort_type eq "-"){
		foreach my $word (sort {$a cmp $b} keys %hash) {
			foreach my $sum (sort {$hash{$word}{$b} <=> $hash{$word}{$a}} keys %{$hash{$word}}) {
				push @{$ref_sorted},$word.$hash{$word}{$sum};
			}
		}
	}
}

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
Program:
Version: $version
Contact:zeng huaping<zenghp\@biomarker.com.cn> 

Usage:			
  Options:
  -id      <str>    dir of inputfile(xxx.map xxx.pwd xxx.loc),forced
  -k       <str>    key of output file（不能以数字开头）,                       forced
  -igeno   <file>   genotype file (.genotype),                 forced
  -od      <str>    dir of output file,                       default keyResult
  -r       <float>  ratio of noise point,                     default 0.05
  --mapEva          mapEvaluation Flag

  -p       <int>    process                                   default 50
  -q       <str>    queue                                     default general.q
  -h         Help

USAGE
	print $usage;
	exit;
}
