#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::RankCorrelation;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fMap,$fAlign,$fKey,$dOut,$PNG,$x_title,$y_title,$width,$height);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$fMap,
				"o:s"=>\$dOut,
				"k:s"=>\$fKey,
				"x:s"=>\$x_title,
				"y:s"=>\$y_title,
				"w:s"=>\$width,
				"h:s"=>\$height,
				"png"=>\$PNG,
				) or &USAGE;
&USAGE unless ($fMap  and $fKey);

#die ` pod2text $0 ` unless ($fMap and $fAlign and $fChroLength and $fKey) ;

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
my $ConvertSVGtoPNG="$Bin/svg2xxx_release/svg2xxx";

$x_title||= "linkage groups";
$y_title||= "genome";
$width||=2000;
$height||=1600;
#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%Map,%Marker,%Alignment,%chroLength) = ();

#-------------------------------------------------------------------
# Get Map info
#-------------------------------------------------------------------

open (IN,"<",$fMap) or die $!;
$/="\n";
my ($group,$order);
my %spearman;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^$/ || /^;/ || /^\#/) ;
	if (/^group\s+(\d+)/) {
		$group=$1;
		$order=0;
	}else{
		my ($marker,$cm,undef)=split(/\s+/,$_);
		$Map{$group}{$marker}{"order"}=$order;
		$Map{$group}{$marker}{"cm"}=$cm;
		$Marker{$marker}{"G"}{"group"}=$group;
		$Marker{$marker}{"G"}{"order"}=$order;
		$Marker{$marker}{"G"}{"cm"}=$cm;
		my ($chr,$pos)=split(/\_/,$marker);
		$chr=~s/\D//g;
		$Marker{$marker}{"C"}{$chr}=$pos;
		$Alignment{$chr}{$marker}=$pos;
		$chroLength{$chr}=$pos if (!exists $chroLength{$chr} || $chroLength{$chr} < $pos);
		$order++;
		$spearman{$group}{$chr}++;
	}
}
close (IN) ;
open Out,">$dOut/$fKey.spearman.xls";
print Out "#LGID\tMarkerNum\tMaxChr\tSpearman\n";
foreach my $chr (sort keys %spearman) {
	my $group=(sort{$spearman{$chr}{$b}<=>$spearman{$chr}{$a}} keys %{$spearman{$chr}})[0];
	my @marker1=sort{$Map{$chr}{$a}<=>$Map{$chr}{$b}}keys %{$Map{$chr}};
	my @marker2=sort{$Alignment{$group}{$a}<=>$Alignment{$group}{$b}} keys %{$Alignment{$group}};
	my (@m1,@m2);
	for (my $i=0;$i<@marker1;$i++) {
		push @m1,$i;
		push @m2,grep{$marker2[$_] eq $marker1[$i]} 0..$#marker2;
	}		
	my $c = Statistics::RankCorrelation->new( \@m1, \@m2,   );
	my $n=$c->spearman;
	print Out $chr,"\t",scalar @marker1,"\t",$chr,"\t",$n,"\n";
}
close Out;

#print Dumper %Map;die;

# ------------------------------------------------------------------
# Get alignment info
# ------------------------------------------------------------------
#my %chroLength;
#print Dumper %Alignment;die;
#------------------------------------------------------------------
# Get chro length
#------------------------------------------------------------------

#print Dumper %chroLength;die;

# ------------------------------------------------------------------
# Draw SVG parameter
# ------------------------------------------------------------------
my $LeftMargin=50;
my $RightMargin=50;
my $TopMargin=50;
my $BottomMargin=50;

my $TextColor="black";
my $TextSize="30";
my $TextLineHeight=20;
my $LineWidth=2;

my $plotAreaWidth = $width;
my $plotAreaHeight = $height;

## color 

my %color = ();
foreach my $chro (keys %chroLength) {

	my $blue = int(rand(255)) + 1;
	my $green = int(rand(255)) + 1;
	my $red = int(rand(255)) + 1;

	$color{$chro} = "rgb($blue,$green,$red)";
}

my $axisColor = "#525252";
my $gridColor = "#D6D6D6";
# ------------------------------------------------------------------
# SVG viewer
# ------------------------------------------------------------------

open (SVG,">$dOut/$fKey.svg") or die $!;

my $yAxis = sum(values %chroLength);
my $xAxis = 0;
my %mapDist = ();

foreach my $lg (keys %Map) {

	my @order = sort {$Map{$lg}{$a}{'order'} <=> $Map{$lg}{$b}{'order'}} keys %{$Map{$lg}};
	$xAxis += $Map{$lg}{$order[-1]}{'cm'} - $Map{$lg}{$order[0]}{'cm'};
	$mapDist{$lg} = $Map{$lg}{$order[-1]}{'cm'} - $Map{$lg}{$order[0]}{'cm'};
}

my $xScale = $xAxis / $plotAreaWidth;
my $yScale = $yAxis / $plotAreaHeight;

my $axisWidth = 10;
my $x_axisMargin = 30;
my $y_axisMargin = 10;

my $PaperHeight = $TopMargin + $BottomMargin + $plotAreaHeight + $LineWidth*(scalar keys %Alignment) + 2*$TextLineHeight + $axisWidth + $x_axisMargin;
my $PaperWidth = $LeftMargin + $RightMargin + $plotAreaWidth + $LineWidth*(scalar keys %Map) + 2*$TextLineHeight + $axisWidth + $y_axisMargin;

my $x = 0;
my $y = 0;

print SVG &svg_paper($PaperWidth,$PaperHeight,"white"),"\n";
print SVG &svg_txt($LeftMargin,($TopMargin+$plotAreaHeight)/2,$TextSize,$TextColor,$y_title,3),"\n";## draw y axis title

$x= $LeftMargin + 2*$TextLineHeight+$y_axisMargin;
$y = $TopMargin;

### draw from the top to the bottom
my @chro=keys %Alignment;
my %stat;
for (my $i=0;$i<@chro;$i++) {
	my $chro_num=$chro[$i];
	$chro_num=~s/\D+//g;
	$stat{$chro[$i]}=$chro_num;
}

foreach my $chro (sort {$stat{$a}<=>$stat{$b}} keys %stat) {

	print SVG &svg_txt($LeftMargin+$TextLineHeight,$y+$chroLength{$chro}/(2*$yScale),$TextSize,$TextColor,$chro),"\n";

	## draw y axis,grid,scatter point of curr_chro

	print SVG &svg_rect($x,$y+$LineWidth,$axisWidth,$chroLength{$chro}/$yScale,$axisColor),"\n";  ## y axis 
	
	print SVG &svg_line($x+$axisWidth,$y,$x+$axisWidth+$plotAreaWidth+$LineWidth*(scalar keys %Map),$y,$gridColor,$LineWidth,0.8),"\n"; ## grid line horizontal
	my %lg = map {$_,1} (map {$Marker{$_}{"G"}{"group"}} keys %{$Alignment{$chro}}); 
	foreach my $lg (sort {$a <=> $b} keys %lg) {
		my @order = sort {$Map{$lg}{$a}{'order'} <=> $Map{$lg}{$b}{'order'}} keys %{$Map{$lg}};
		my @marker = @{&intersection([keys %{$Alignment{$chro}}],\@order)};
		#Linkage sort by Alighment
		my $flag1=0;
		my $flag2=0;
		for (my $i=0;$i<@order-1;$i++) {
			next if (!exists $Alignment{$chro}{$order[$i]}||!exists $Alignment{$chro}{$order[$i+1]});
			if ($Alignment{$chro}{$order[$i]} < $Alignment{$chro}{$order[$i+1]}){
				$flag1++;
			}else{
				$flag2++;
			}
		}
		my $radius = 2;

		my $grid_x = $x+$axisWidth + sum(map {$mapDist{$_}} grep {$_ < $lg} keys %Map)/$xScale + $LineWidth *(scalar grep {$_ < $lg} keys %Map);
		my $grid_y = $y;
		foreach my $marker (@marker) {
			if ($flag2 > $flag1) {
				$x = $grid_x + ($Map{$lg}{$order[$#order]}{'cm'} - $Map{$lg}{$marker}{'cm'})/$xScale;
			}else{
				$x = $grid_x + ($Map{$lg}{$marker}{'cm'} - $Map{$lg}{$order[0]}{'cm'})/$xScale;
			}
			$y = $grid_y + $Alignment{$chro}{$marker} / $yScale ;
			if (!exists $color{$chro}) {
				die;
			}
			print SVG &svg_circle($x,$y,$radius,$color{$chro}),"\n";

			$x = $LeftMargin + 2*$TextLineHeight+$y_axisMargin;
			$y = $grid_y;

		}
	}
	
	$x= $LeftMargin + 2*$TextLineHeight+$y_axisMargin;
	$y += $chroLength{$chro}/$yScale + $LineWidth;

}

#$y = $TopMargin + $LineWidth*(scalar keys %Alignment) + $plotAreaHeight;

#$y = $plotAreaHeight;
$x = $LeftMargin + 2*$TextLineHeight + $y_axisMargin + $axisWidth ;


foreach my $lg (sort {$a <=> $b} keys %Map) {

	print SVG &svg_txt(($LeftMargin+$plotAreaWidth)/2,$y + $axisWidth + 2*$x_axisMargin+$TextLineHeight,$TextSize,$TextColor,$x_title),"\n";

	print SVG &svg_txt($x+$mapDist{$lg}/(2*$xScale),$y+$axisWidth+$x_axisMargin,$TextSize,$TextColor,$lg),"\n"; ## x axis label
	
	print SVG &svg_rect($x,$y,$mapDist{$lg}/$xScale,$axisWidth,$axisColor),"\n";  ## x axis 

	print SVG &svg_line($x+$mapDist{$lg}/$xScale,$y,$x+$mapDist{$lg}/$xScale,$TopMargin,$gridColor,$LineWidth,0.8),"\n"; ## grid line vertical

	$x += $mapDist{$lg}/$xScale + $LineWidth;
#	$y = $TopMargin + $LineWidth*(scalar keys %Alignment) + $plotAreaHeight;
}

print SVG &svg_end(),"\n";

close (SVG) ;
`cairosvg $dOut/$fKey.svg -o $dOut/$fKey.png --dpi 300`;
`cairosvg $dOut/$fKey.svg -o $dOut/$fKey.pdf `;
## convert to png 
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub svg_paper (){#&svg_paper(width,height,[color])
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
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

sub svg_rect () {#&svg_rect(x,y,width,height,color,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
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
	my $line="<circle cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" r=\"$svg_x[2]\" stroke=\"$svg_x[3]\" stroke-width=\"0\" fill=\"$svg_x[3]\"/>";
	return $line;
}

sub svg_path () {#colorfill,colorstroke,strokewidth,coloropacity,$path
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $width=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_path=shift(@svg_x);
	my $line="<path d= \"$svg_path\" fill=\"$svg_color\" stroke=\"$svg_color2\" stroke-width=\"$width\" opacity=\"$svg_trans\"/>\n";
	return $line;
}

sub chroLen {#
	my ($file,$ref)=@_;
	open (IN,"<",$file) or die $!;
	my ($chro,$sum,$count);
	while (<IN>) {
		next if (/^$/) ;
		if (/^>/) {
			if ($chro) {
				$sum+=$count;
				$ref->{$chro}=$count;
			}
			($chro)=($_=~/^>(\w+)/);
			$count=0;
			next;
		}else{
			s/\s//g;
			$count+=length($_);
		}
	}
	$sum+=$count;
	$ref->{$chro}=$count;
	close (IN) ;
	return $sum;
}
sub max{#&max(lists or arry);
	#���б��е����ֵ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

sub intersection {#
	my ($A,$B)=@_;
	my %uniqA=map {$_,1} @{$A};
	my %uniqB=map {$_,1} @{$B};
	my %merge=();
	my %overlap=();
	foreach  (keys %uniqA,keys %uniqB) {
		$merge{$_}++ && $overlap{$_}++;
	}
	my @result = keys %overlap;
	return \@result;
}

sub sum {#
	my $sum=0;
	$sum+=$_  foreach (@_) ;
	return $sum;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>

	map file: 
	  group X
	  MarkerID	MapDistance

	alignment file:
	  MarkerID	ChromosomeID	startPos of alignment

	info file:
	  ChromosomeID	ChromosomeLength

Usage:
  Options:
  -m <file>  map file, JoinMap map format, forced.
  -x <str>	 title of x axis ,default [linkage groups], optional.
  -y <str>   titile of y axis ,default [genome], optional.
  -w <str>   width of svg paper ,default [2500], optional.
  -h <str>   height of svg paper ,default [1600], optional.
  -k <str>   key of output file , forced.
  -o <dir>   direction where restult proced.
  -png       Convert svg to png.

USAGE
	print $usage;
	exit;
}
