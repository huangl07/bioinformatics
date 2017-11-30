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
my ($fIn1,$fIn2,$fIn3,$outdir,$fOut,$fOut1,$Reverse,$switch);
GetOptions(
				"help|?" =>\&USAGE,
				"i1:s"=>\$fIn1,
				"i2:s"=>\$fIn2,
				"i3:s"=>\$fIn3,
				"o:s"=>\$fOut,
				"o1:s"=>\$fOut1,
				"d:s"=>\$outdir,
				"r:s"=>\$Reverse,
				"s:s"=>\$switch,
				) or &USAGE;
$switch||="no";
&USAGE unless ($switch eq "run" and $fIn1 and $fIn2 and $fIn3);

#------------------------------------------------read in files


#print Dumper %all;die;
$Reverse = 0 unless (defined $Reverse) ;
$outdir||="./";
`mkdir $outdir` unless (-d $outdir); 
my $workdir=`pwd`;chomp $workdir;

#------------------------------------------------------------------
#  global  variable
#------------------------------------------------------------------
my %Maps;

#------------------------------------------------------------------
# Get map info
#------------------------------------------------------------------

get_map($fIn1, \%{$Maps{'1'}{'left'}});
get_map($fIn2, \%{$Maps{'1'}{'middle'}});
get_map($fIn3, \%{$Maps{'1'}{'right'}});

#print Dumper %Maps;die;

# ------------------------------------------------------------------
# Draw SVG parameter
# ------------------------------------------------------------------
my $LeftMargin=50;
my $RightMargin=50;
my $TopMargin=50;
my $BottomMargin=50;

my $TextColor="black";
my $TextSize="12";
my $TextWidth=20;
my $TextSizeX="14";

my $TitleHeight=10;
my $titleHeight=10;
my $Title_title_space=20;
my $title_chro_space=10;

my $MarkerTxtSpace=45;
my $MarkerLineSpace=20;
my $Marker_LineSpace=40;
my $CmTxtSpace=60;
#my $CmLineSpace=30;

my $LineSize=1;
my $LineColor="#2BD5B3";
#my $LineColor="#6B8E23";
my $ShortLineColor="black";

my $space_betw_chro=120;

my $ChroWidth=20;
my $ChroColor="#999900";
my $factor=10;
my $MiddleColor="#EE7942";
my $LeftColor="#3CB371";
my $RightColor="#1C86EE";
my $filterColor = "#FFFFFF"; 

#-------------------------------------------------得到三个文件中相应group的最大遗传距离，并存到一个数组中    
my %max_chro_length;
foreach my $group (sort keys %Maps) { 
	my $max=0;
	foreach my $mapi (qw(left middle right)) { 
		my ($end) = sort {$Maps{$group}{$mapi}{$b}{'order'} <=> $Maps{$group}{$mapi}{$a}{'order'}} keys %{$Maps{$group}{$mapi}};
		 $max=$Maps{$group}{$mapi}{$end}{'cm'} ;
		 push @{$max_chro_length{$group}},$max;
	}
	unshift  @{$max_chro_length{$group}},"0" if (!defined $Maps{$group}{'left'});
	push   @{$max_chro_length{$group}},"0" if (!defined $Maps{$group}{'right'});
} 

#print Dumper %max_chro_length;
#die;

my ($pos,$synteny)=&get_colinearity_data(\%Maps);

my %pos=%$pos;my %synteny=%$synteny;
foreach my $groups (sort keys %Maps) {
	print STDOUT "LG$groups is drawing......\n";
	my $png = defined $fOut ? "$fOut.$groups":$groups;  
	open (SVG,">","$outdir/$png.FAM.svg") or die $!;
	my @tem=@{$max_chro_length{$groups}};
	my @temp = sort {$a <=> $b} @tem;
#	print @tem,"\n";die;
	my $max_length=$temp[-1];
	$factor= (keys %{$Maps{$groups}{'middle'}})*$TextSizeX >=$max_length*$factor ? ((keys %{$Maps{$groups}{'middle'}})*$TextSizeX)/$max_length:$factor;
	$factor= (keys %{$Maps{$groups}{'left'}})*$TextSizeX >=$tem[0]*$factor ? ((keys %{$Maps{$groups}{'left'}})*$TextSizeX)/$tem[0]:$factor;
	$factor= (keys %{$Maps{$groups}{'right'}})*$TextSizeX >=$tem[2]*$factor ? ((keys %{$Maps{$groups}{'right'}})*$TextSizeX)/$tem[2]:$factor;
	
	###  draw SVG paper 
	my $PaperWidth=$LeftMargin+$RightMargin+($MarkerTxtSpace*3+$MarkerLineSpace*6+$ChroWidth*3+$space_betw_chro*2)+$CmTxtSpace*3+$Marker_LineSpace*6;
	my $PaperHeight=$TopMargin+$BottomMargin+$max_length*$factor+$ChroWidth+$TitleHeight+$titleHeight+$Title_title_space+$title_chro_space;
	print SVG &svg_paper($PaperWidth,$PaperHeight);
	print SVG "<defs>\n<linearGradient  id=\"left_chro_bar_filter\"  x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$LeftColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$filterColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$LeftColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";
	print SVG "<defs>\n<linearGradient  id=\"middle_chro_bar_filter\"  x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$MiddleColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$filterColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$MiddleColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";
	print SVG "<defs>\n<linearGradient  id=\"right_chro_bar_filter\"  x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$RightColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$filterColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$RightColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";
	print SVG &svg_txt($PaperWidth/2,$TopMargin+$TitleHeight,$TextSize,$TextColor,"LG $png");
	my $x1=$LeftMargin+$MarkerLineSpace+$MarkerTxtSpace+$Marker_LineSpace;
	my $x2=$LeftMargin+$MarkerLineSpace*3+$MarkerTxtSpace*2+$space_betw_chro+$ChroWidth+$CmTxtSpace+$Marker_LineSpace*3;
	my $x3=$LeftMargin+$MarkerLineSpace*5+$MarkerTxtSpace*3+$space_betw_chro*2+$ChroWidth*2+$CmTxtSpace*2+$Marker_LineSpace*5;

	### draw chro bar
	print SVG &svg_txt($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space + ($max_length - $tem[0]) * $factor / 2,$TextSize,$TextColor,"$png.female");
	print SVG &svg_circle($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space+($max_length - $tem[0]) * $factor / 2,$ChroWidth/2,$LeftColor,"left_chro_bar_filter"),"\n";
	print SVG &svg_circle($x1+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[0])*$factor+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2, $ChroWidth/2, $LeftColor,"left_chro_bar_filter"),"\n"; 
	print SVG &svg_rect($x1,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$ChroWidth,$tem[0]*$factor,$LeftColor,"left_chro_bar_filter") ;
	
	print SVG &svg_txt($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space + ($max_length - $tem[1]) * $factor/2,$TextSize,"$TextColor","$png.sexAver");
	print SVG &svg_circle($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1])* $factor/2,$ChroWidth/2,$MiddleColor,"middle_chro_bar_filter"),"\n";
	print SVG &svg_circle($x2+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[1])*$factor+$Title_title_space+$title_chro_space + ($max_length - $tem[1])* $factor/2,$ChroWidth/2,$MiddleColor,"middle_chro_bar_filter"),"\n"; 
	print SVG &svg_rect($x2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1])* $factor/2,$ChroWidth,$tem[1]*$factor,$MiddleColor,"middle_chro_bar_filter") ;
	
	print SVG &svg_txt($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$Title_title_space + ($max_length - $tem[2])* $factor/2,$TextSize,"$TextColor","$png.male");
	print SVG &svg_circle($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2])* $factor/2,$ChroWidth/2,$RightColor,"right_chro_bar_filter"),"\n";
	print SVG &svg_circle($x3+$ChroWidth/2,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+($tem[2])*$factor+$Title_title_space+$title_chro_space + ($max_length - $tem[2])* $factor/2,$ChroWidth/2,$RightColor,"right_chro_bar_filter"),"\n"; 
	print SVG &svg_rect($x3,$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2])* $factor/2,$ChroWidth,$tem[2]*$factor,$RightColor,"right_chro_bar_filter") ;
	
	### draw locus info and map position info 
	foreach my $mapi (qw(left middle right)) {
	
		my @Marker=keys %{$Maps{$groups}{$mapi}};
		my $k=0;
		foreach my $marker(sort {$Maps{$groups}{$mapi}{$a}{'order'} <=> $Maps{$groups}{$mapi}{$b}{'order'}} keys  %{$Maps{$groups}{$mapi}} ){
			my $y=$Maps{$groups}{$mapi}{$marker}{'cm'}*$factor;
			my $Y=$Maps{$groups}{$mapi}{$marker}{'cm'};

			if ($mapi eq 'left') {
				my $length=$tem[0]/(scalar @Marker)*$factor;
				my $MY=$length*$k;
				## locus
				print SVG &svg_txt($LeftMargin,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$TextSize,"$TextColor","$marker");
				print SVG &svg_line($x1-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x1,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$ShortLineColor,$LineSize);
				print SVG &svg_line($x1-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x1-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$ShortLineColor,$LineSize);
				## map position 
				print SVG &svg_txt($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$TextSize,"$TextColor","$Y cm");
				print SVG &svg_line($x1+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x1+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$ShortLineColor,$LineSize);
				print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$ShortLineColor,$LineSize);
				$pos{$groups}{'left'}{$marker}=$MY;
				$k++;
			}elsif($mapi eq 'middle') {
				my $length;
				$length=$tem[1]/(scalar @Marker)*$factor;
				my $MY=$length*$k;
				if ($pos{$groups}{"7k"} eq "reverse") {
					my $bottomY=$tem[1]*$factor+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2;
					print SVG &svg_txt($x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x2-$MarkerLineSpace,$bottomY-$y,$x2,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2-$MarkerLineSpace-$Marker_LineSpace,$bottomY-$MY,$x2-$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$Y cm");
					print SVG &svg_line($x2+$ChroWidth,$bottomY-$y,$x2+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$ShortLineColor,$LineSize);
					$pos{$groups}{"reverse_middle"}{$marker}=$bottomY-$MY;
				}else{
					print SVG &svg_txt($x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x2-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x2,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x2-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$TextSize,"$TextColor","$Y cm");
					print SVG &svg_line($x2+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x2+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$ShortLineColor,$LineSize);
					print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$ShortLineColor,$LineSize);
					$pos{$groups}{"non_reverse_middle"}{$marker}=$MY;
				}
				$k++;
			}elsif ($mapi eq 'right'){
				my $length;my $chrolength;
				$length=$tem[2]/(scalar @Marker)*$factor;
				$chrolength=$tem[2];
				my $MY=$length*$k;
				if ($pos{$groups}{"9k"} eq "reverse") {
					my $bottomY=$chrolength*$factor+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2;
					print SVG &svg_txt($x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x3-$MarkerLineSpace,$bottomY-$y,$x3,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3-$MarkerLineSpace-$Marker_LineSpace,$bottomY-$MY,$x3-$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$TextSize,"$TextColor","$Y cm");
					print SVG &svg_line($x3+$ChroWidth,$bottomY-$y,$x3+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3+$ChroWidth+$MarkerLineSpace,$bottomY-$y,$x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$bottomY-$MY,$ShortLineColor,$LineSize);
					$pos{$groups}{"reverse_right"}{$marker}=$bottomY-$MY;
				}else{
					print SVG &svg_txt($x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$TextSize,"$TextColor","$marker");
					print SVG &svg_line($x3-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$x3,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3-$MarkerLineSpace-$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2, $x3-$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$ShortLineColor,$LineSize);

					print SVG &svg_txt($x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$TextSize,"$TextColor","$Y cm");
					print SVG &svg_line($x3+$ChroWidth,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$x3+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$ShortLineColor,$LineSize);
					print SVG &svg_line($x3+$ChroWidth+$MarkerLineSpace,$y+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$x3+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace,$MY+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$ShortLineColor,$LineSize);
					$pos{$groups}{"non_reverse_right"}{$marker}=$MY;
				}
				$k++;
			}
		}
	}
	### draw synteny line  
	my @marker1=keys %{$Maps{$groups}{'left'}};
	my @marker2=keys %{$Maps{$groups}{'middle'}};
	foreach my $marker1 (@marker1) {
		if (defined $pos{$groups}{"reverse_middle"}{$marker1}) {
			print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{'left'}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"reverse_middle"}{$marker1},$LineColor,$LineSize);
		}elsif (defined $Maps{$groups}{'middle'}{$marker1}) {
			print SVG &svg_line($x1+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{'left'}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[0]) * $factor / 2,$x2-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"non_reverse_middle"}{$marker1}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2 ,$LineColor,$LineSize);
		}
	}
	if ($pos{$groups}{"7k"} eq "reverse") {
		foreach my $marker2 (@marker2) {
			if (defined $pos{$groups}{"reverse_right"}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"reverse_middle"}{$marker2},$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"reverse_right"}{$marker2},$LineColor,$LineSize);
			}elsif (defined $Maps{$groups}{'right'}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"reverse_middle"}{$marker2},$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"non_reverse_right"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$LineColor,$LineSize);
			}
		}
	}else{
		foreach my $marker2 (@marker2) {
			if (defined $pos{$groups}{"reverse_right"}{$marker2}) {
				print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"non_reverse_middle"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"reverse_right"}{$marker2},$LineColor,$LineSize);
			}elsif (defined $Maps{$groups}{'right'}{$marker2}) {
				 print SVG &svg_line($x2+$ChroWidth+$MarkerLineSpace+$Marker_LineSpace+$CmTxtSpace,$pos{$groups}{"non_reverse_middle"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[1]) * $factor / 2,$x3-$MarkerLineSpace-$MarkerTxtSpace-$Marker_LineSpace-4,$pos{$groups}{"non_reverse_right"}{$marker2}+$TopMargin+$TitleHeight+$titleHeight+$ChroWidth/2+$Title_title_space+$title_chro_space + ($max_length - $tem[2]) * $factor / 2,$LineColor,$LineSize);
			}
		}
	}
	print SVG &svg_end ();
	close (SVG);
	chdir $outdir;
	`perl $Bin/svg2xxx_release/svg2xxx  "$png.FAM.svg" -t -w $PaperWidth -h $PaperHeight `;
	chdir $workdir;
}
foreach my $groups (sort keys %synteny) {
	my $png = defined $fOut ? $fOut:$groups;  
	open (STAT,">$outdir/$png.FAM.synteny.stat.xls") or die $!;
	print STAT "LG$png\n";
	print STAT"${$synteny{$groups}}[0]\/${$synteny{$groups}}[2]\t${$synteny{$groups}}[1]\/${$synteny{$groups}}[3]\n";
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub get_map {#
	my ($fIn, $ref_map) = @_;
	my $order = 0;
	open (IN,$fIn) or die $!;
	$/="\n";
	while (<IN>) {
		chomp;
		s/\r//g;
		next if /^\s/ || /^;/ || /^group/;
		next unless (/^(\S+)\s+(\S+)/) ;
		$ref_map->{$1}{'order'} = $order++ ;
		$ref_map->{$1}{'cm'} = $2 ;
	} 
	close(IN);
	## check the second map file
	die "The second map file $fIn is empty, please check" if (scalar keys %{$ref_map} == 0) ;
}

sub get_colinearity_data {#  共线性数量，公共marker数量的计算及染色体是否要翻转的判断
	my %copyMaps=%{$_[0]};
	my %position;
	my %synteny_stat;
	foreach my $groups (sort keys %copyMaps) {
		if ($Reverse==0) {
			my ($synteny_stat1,$synteny_stat2,$a1,$a2);
			my (@order1,@order11,@order2,@order22);
			my $k=0;my $i=0;my $j=0;my $p=0;
			my @marker2=(sort {$copyMaps{$groups}{'middle'}{$a}{'order'} <=> $copyMaps{$groups}{'middle'}{$b}{'order'}} keys  %{$copyMaps{$groups}{'middle'}} );
			my @Rmarker2=reverse @marker2;
			my @marker1=(sort {$copyMaps{$groups}{'left'}{$a}{'order'} <=> $copyMaps{$groups}{'left'}{$b}{'order'}} keys  %{$copyMaps{$groups}{'left'}} );
			my @marker3=(sort {$copyMaps{$groups}{'right'}{$a}{'order'} <=> $copyMaps{$groups}{'right'}{$b}{'order'}} keys  %{$copyMaps{$groups}{'right'}} );
			my @Rmarker3=reverse @marker3;
			foreach my $marker2 (@marker2) {
				 $k++;
				$position{$groups}{'middle'}{$marker2}=$k;
				
			}
			foreach my $Rmarker2 (@Rmarker2) {
				 $j++;
				$position{$groups}{"77"}{$Rmarker2}=$j;
				
			}
			foreach my $marker3 (@marker3) {
				 $p++; 
				$position{$groups}{'right'}{$marker3}=$p;
			
			}
			foreach my $Rmarker3 (@Rmarker3) {
				 $i++;
				$position{$groups}{"99"}{$Rmarker3}=$i;
			
			}
			foreach my $marker1 (@marker1) {
				if (defined $position{$groups}{'middle'}{$marker1}) {
					push  @order1,$position{$groups}{'middle'}{$marker1};
				}
				if (defined $position{$groups}{"77"}{$marker1}) {
					push @order11,$position{$groups}{"77"}{$marker1};
				}
			}
			$a1=@order1;
			if (scalar @order1==0 || scalar @order11==0){
				$position{$groups}{"7k"}="not reverse";
			}else{
				my @result1=PickMaxGroup(@order1);
				my @result11=PickMaxGroup(@order11);

		#print scalar @result1,"\n",scalar @result11,"\n";
				if (scalar @result1 <  scalar @result11) {
					$synteny_stat1=scalar @result11;
					$position{$groups}{"7k"}="reverse";
					@marker2=@Rmarker2;
				}else{
					$position{$groups}{"7k"}="not reverse";
					$synteny_stat1=scalar @result1;
				}	
			}
			foreach my $marker2 (@marker2) {
				if (defined $position{$groups}{'right'}{$marker2}) {
				push  @order2,$position{$groups}{'right'}{$marker2};
				}
				if (defined $position{$groups}{"99"}{$marker2}) {
				push  @order22,$position{$groups}{"99"}{$marker2};
				}
			}
			$a2=@order2;
			if (scalar @order2==0 || scalar @order22==0 ) {
				$position{$groups}{"9k"}="not reverse";
			}else{
				my  @result2=PickMaxGroup(@order2);
				my  @result22=PickMaxGroup(@order22);
				if (scalar @result2 <  scalar @result22) {
					$synteny_stat2=scalar @result22;
					$position{$groups}{"9k"}="reverse";
				}else{
					$position{$groups}{"9k"}="not reverse";
					$synteny_stat2=scalar @result2;
				}	
			}
			$synteny_stat1= !defined $synteny_stat1 ? 0:$synteny_stat1; 
			$synteny_stat2= !defined $synteny_stat2 ? 0:$synteny_stat2;
			$a1= !defined $a1 ? 0:$a1;
			$a2= !defined $a2 ? 0:$a2;
			push @{$synteny_stat{$groups}},$synteny_stat1;
			push @{$synteny_stat{$groups}},$synteny_stat2;
			push @{$synteny_stat{$groups}},$a1;
			push @{$synteny_stat{$groups}},$a2;
		}else{
			$position{$groups}{"9k"}="not reverse";
			$position{$groups}{"7k"}="not reverse";
		}
		
	}
	return (\%position,\%synteny_stat);
}

sub PickMaxGroup
{
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

sub svg_paper (){#&svg_paper(width,height,[color])
	my $svg_drawer = "litc"."@"."biomarker\.com\.cn";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub svg_end (){#
	return "</svg>\n";
}

sub svg_txt (){#&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_line (){#&svg_line(x1,y1,x2,y2,color,width,[opacity])
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" opacity=\"$svg_x[6]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_polyline (){#colorfill,colorstroke,width,\@point
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_width=shift(@svg_x);
	my $svg_points=join(" ",@{$svg_x[-1]});
	my $line="<polyline fill=\"$svg_color\" stroke=\"$svg_color2\" stroke-width=\"$svg_width\" points=\"$svg_points\"/>\n";

	#print "$line\n";
	return $line;

	#<polyline points="0,0 0,20 20,20 20,40 40,40 40,60" style="fill:white;stroke:red;stroke-width:2"/>
}

sub svg_rect () {#&svg_rect(x,y,width,height,color,Id,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=1;
	}
	my $line="<rect style=\"fill:url(#$svg_x[5])\" x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[6]\"/>\n";
	return $line;
}

sub svg_rect1 () {#&svg_rect(x,y,width,height,color,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[6]) {
		$svg_x[6]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[6]\"/>\n";
	return $line;
}
sub svg_polygon () {#colorfill,colorstroke,coloropacity,point1,point2,...
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_ellipse () {#&svg_ellipse(cx,cy,rx,ry,colorfill,colorstroke,width,[coloropacity])
	my @svg_x=@_;
	my $line= "<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\"/>\n";
	if (defined $svg_x[7]) {
		$line="<ellipse cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" rx=\"$svg_x[2]\" ry=\"$svg_x[3]\" fill=\"$svg_x[4]\" stroke=\"$svg_x[5]\" stroke-width=\"$svg_x[6]\" opacity=\"$svg_x[7]\"/>\n";
	}
	return $line;
}

sub svg_circle () {#&svg_circle(cx,cy,r,color)
	my @svg_x=@_;
	my $line="<circle style=\"fill:url(#$svg_x[4])\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" r=\"$svg_x[2]\" stroke=\"$svg_x[3]\" stroke-width=\"0\" fill=\"$svg_x[3]\"/>";
	return $line;
}


sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:draw_ThreeGroup.pl
Version: version1.1
Contact: Xin HuaiGen <xinhg\@biomarker.com.cn>
Description:you can use it like:
	        perl draw_ThreeGroup.pl  -s run
	        or:
	        perl draw_ThreeGroup.pl  -i1 *.female.map -i2 *.sexAver.map -i3 *.male.map  -o xxx  -s run

Usage:map.pl <可用于画三个图谱间的共线性关系，从左往右依次是雌中雄,运行时不输入文件，会在当前路径下找形如*.female.map、*.male.map、*.sexAver.map的文件，并绘图>
  Options:
  -i1 <file>  *.map file1,format,with less markers,usuMapsy the *.female.map
  -i2 <file>  *.map file2,format,whth most markers,usuMapsy the *.sexAver.map
  -i3 <file>  *.map file3,format,whit less markers,usuMapsy the *.male.map
  -o <string>  key of output file,if not given,output file would be named after group number 
  -d <string>  output dir, default,"synteny"
  -r <number>  optional,0 means judge whether reverse,any other numbers mean don\'t reverse straightly ,default,0
  -s <control> must be given as 'run' or it won't run
  -h         Help

USAGE
	print $usage;
	exit;
}														
