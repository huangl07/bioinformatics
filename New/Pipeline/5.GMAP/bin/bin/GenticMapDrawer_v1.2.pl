#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.5.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($Map,$Svg,$Row,$Column,$PNG);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$Map,
				"o:s"=>\$Svg,
				"r:s"=>\$Row,
				"c:s"=>\$Column,
				) or &USAGE;
&USAGE unless ($Map and $Svg and $Column and $Row);



# ------------------------------------------------------------------
# Get Map info
# ------------------------------------------------------------------
my %SLAF;  my @files;
my %MARKER;
my $Chro;
if (-f $Map) {
	push @files,$Map;
}elsif(-d $Map){
	@files=glob("$Map/*.map");
}



foreach my $file (@files) {
	open (IN,"<",$file) or die $!;
	if ($file=~/(\D)\.map/) {
		($Chro)=$file=~/(\d+)(\.|\_)(\w+)\.map/;
		#print "$Chro\n";
	}
	elsif ($file=~/(\d+).map/) {
		($Chro)=$file=~/(\d+)\.map/;

	}
	else{
		print "The mapfile must be named by order , please check your mapname.";
		die;
	}


	while (<IN>){
		chomp;
		s/\r//g;
		next if (/^$/ || /^;/ || /^\#/ || /^nloc/) ;
		if (/^group\s+(\w+)/i) {
			$Chro=$1;
		}else{
			my ($Id,$Start)=$_=~/^(\S+)\s+(\S+)/;
			$MARKER{$Chro}{$Id}=$_;
			
		}
	}
	foreach my $key (sort keys %{$MARKER{$Chro}}) {
		my ($Id,$Start)=$MARKER{$Chro}{$key}=~/^(\S+)\s+(\S+)/;
		$SLAF{$Chro}{$Id}{"end"}=$Start;
		$SLAF{$Chro}{$Id}{"id"}=$Id;
#		print "$Start\n";
	}


	close (IN) ;
}




# print Dumper %SLAF;die; 
# ------------------------------------------------------------------
# Get len of reference
# ------------------------------------------------------------------
my %Len;

foreach my $chro (keys %SLAF) {
	my $len=(sort {$SLAF{$chro}{$b}{"end"} <=> $SLAF{$chro}{$a}{"end"}} keys %{$SLAF{$chro}})[0];
	$Len{$chro}=$SLAF{$chro}{$len}{"end"};
}

my ($MAXLEN)=sort {$b <=> $a} values %Len;
#print $MAXLEN;
#die;


# ------------------------------------------------------------------
# Draw SVG parameter
# ------------------------------------------------------------------
my $LeftMargin=20;##左边框空间
my $RightMargin=20;##右边框空间
my $TopMargin=20;##上部空间
my $TopTxtMargin=50;##标头空间
#my $TopTxtMargin=30*$Row;
my $BottomMargin=20;##下部空间

my $WidthSpace=600;##图布宽度空间
my $HeightSpace=800;##图布长度空间

# ------------------------------------------------------------------
my $OneWidthSpace=$WidthSpace/$Column;##单个图宽度空间
my $OneHeightSpace=$HeightSpace/$Row;##单个图长度空间
my $ChroWidth=$OneWidthSpace/11;##染色体宽度或斜线宽度
my $ChroTxtSpace=$ChroWidth/2;##横线长
my $ChroMarkerSpace=$ChroWidth;##Marker宽度
my $ChroCMSpace=$ChroTxtSpace;##遗传距离宽度
my $ScaleCrossLine=50;##线长空间


########my $ChroHeight=$OneHeightSpace/2000*$Markernum;##染色体长度
#
#my $OffsetScale=3;
#my $ScaleHeight=9;##Marker宽
#-------------------------------------------------------------------
#my $ChroTxtSpace=20;##染色体表头空间
#my $ChroWidth=35;##染色体宽
my $ChroColor="#999900";##染色体颜色
#my $ChroInterval=15;##间隔
my $MiddleColor="rgb(247,251,209)";##染色体中间颜色




# ------------------------------------------------------------------
# Drawing
# ------------------------------------------------------------------
my $MaxMarker;
my @Totalnum;
my @Chrnum;
my $TotalChr;
foreach my $chro (sort {$a <=> $b} keys %SLAF) {
	push @Totalnum,(scalar keys %{$SLAF{$chro}});
	push @Chrnum,$chro;
}
$MaxMarker=&max(@Totalnum);
$TotalChr=&max(@Chrnum);
#print "$TotalChr\n\n";die

my @SVGFiles=();
open (SVG,">","$Svg.all.svg") or die $!;
my $PaperWidth=$LeftMargin+$RightMargin+$WidthSpace;#纸张宽度
my $PaperHeight=$TopMargin+$TopTxtMargin*$Row+$BottomMargin+$HeightSpace;#纸张长度


push @SVGFiles,"$Svg.all.svg",$PaperWidth,$PaperHeight;
print SVG &svg_paper($PaperWidth,$PaperHeight),"\n";
print SVG "<defs>\n<linearGradient id=\"Gene\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$MiddleColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";


####################画图

my $n=1;
my $i;my $j;
for ($j=1;$j<= $Row;$j++) {
	for ($i=1;$i<=$Column ;$i++) {
		if ($n<=$TotalChr) {
#			next if (exists %{$SLAF{$n});
				
			my $TotalMarker=(scalar keys %{$SLAF{$n}});##单个连锁群Marker总数
			#print "$TotalMarker\t$OneHeightSpace\n";die;
			my $ChroHeight=($OneHeightSpace*$TotalMarker)/($MaxMarker+50);##染色体长度
			my $Scale=$Len{$n}/($ChroHeight);
			my $ScaleHeight=($ChroHeight/$TotalMarker);
			my $map_title_prefix=(length $n )==1 ?"LG0":"LG";

			print SVG &svg_txt($LeftMargin+$OneWidthSpace/2+$OneWidthSpace*($i-1)-1.4*$ChroWidth,$TopMargin-20*$ScaleHeight+$TopTxtMargin+($TopTxtMargin+$OneHeightSpace)*($j-1),$ChroWidth,"black","$map_title_prefix".$n),"\n";
			print SVG &svg_circle($LeftMargin+$OneWidthSpace/2+$OneWidthSpace*($i-1),$TopMargin+$TopTxtMargin+$ChroWidth/2+($TopTxtMargin+$OneHeightSpace)*($j-1),$ChroWidth/2,$ChroColor),"\n";
			print SVG &svg_circle($LeftMargin+$OneWidthSpace/2+$OneWidthSpace*($i-1),$TopMargin+$TopTxtMargin+$ChroWidth/2+$Len{$n}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),$ChroWidth/2,$ChroColor),"\n";
			print SVG &svg_rect($LeftMargin+$OneWidthSpace/2+$OneWidthSpace*($i-1)-$ChroWidth/2,$TopMargin+$TopTxtMargin+$ChroWidth/2+($TopTxtMargin+$OneHeightSpace)*($j-1),$ChroWidth,$Len{$n}/$Scale,$ChroColor),"\n";

			my $numSLAF=0;
			foreach my $start (sort {$SLAF{$n}{$a}{"end"} <=> $SLAF{$n}{$b}{"end"}} keys %{$SLAF{$n}}) {
				my $posi=$LeftMargin;
				print SVG &svg_line($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)-$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),$posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)-$ChroWidth/2,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),"black",0.03*$ScaleHeight);
				print SVG &svg_line($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)+$ChroWidth/2,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),$posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)+$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),"black",0.03*$ScaleHeight);
				print SVG &svg_txt($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)-3*$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+($numSLAF+0.5)*$ScaleHeight+($TopTxtMargin+$OneHeightSpace)*($j-1),$ScaleHeight,"black",$SLAF{$n}{$start}{"end"}."cM");
				print SVG &svg_txt($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)+2.5*$ChroWidth+3*$ScaleHeight,$TopMargin+$TopTxtMargin+$ChroWidth/2+($numSLAF+0.5)*$ScaleHeight+($TopTxtMargin+$OneHeightSpace)*($j-1),$ScaleHeight,"black",$SLAF{$n}{$start}{"id"});
				print SVG &svg_line($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)-2.5*$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+($numSLAF+0.8)*$ScaleHeight+($TopTxtMargin+$OneHeightSpace)*($j-1),$posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)-$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),"black",0.03*$ScaleHeight);
				print SVG &svg_line($posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)+$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+$SLAF{$n}{$start}{"end"}/$Scale+($TopTxtMargin+$OneHeightSpace)*($j-1),$posi+$OneWidthSpace/2+$OneWidthSpace*($i-1)+2.5*$ChroWidth,$TopMargin+$TopTxtMargin+$ChroWidth/2+($numSLAF+0.8)*$ScaleHeight+($TopTxtMargin+$OneHeightSpace)*($j-1),"black",0.03*$ScaleHeight);
				$numSLAF++;
			}
		}
		else{
			last;
		}
		$n++;
	}
}



	print SVG &svg_end();

close (SVG) ;
`cairosvg $Svg.all.svg -o $Svg.all.png --dpi 300`;
`cairosvg $Svg.all.svg -o $Svg.all.pdf `;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub svg_rect () {#&svg_rect(x,y,width,height,color,[opacity])
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect style=\"fill:url(#Gene)\" x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
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
	my $line="<circle style=\"fill:url(#Gene)\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" r=\"$svg_x[2]\" stroke=\"$svg_x[3]\" stroke-width=\"0\" fill=\"$svg_x[3]\"/>";
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
sub Readfasta {
	my $infile=shift;
	my $seqId_seq=shift;

	my $c=0;
	open IN, $infile || die $!;
	my $seqId;
	while (<IN>) {
		if (/>(\S+)/) {
			$seqId=$1;
			$c++;
		}
		else {
			$_=~s/\s//g;
			$seqId_seq->{$seqId}.=$_;
		}
	}
	close IN;
	return $c;
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
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: SLAFonRef.pl (遗传图谱绘图程序,连锁群小于10个且每个连锁群平均的上图标记<=200时使用)
Version: $version
Contact: Su Yalei <suyl\@biomarker.com.cn> <741311962\@qq.com>
Description：
	Support JoinMap‘s map format files;
	Support Converting svg to png files;
	Every linkage group will be drawn in each file;
	Auto OffSet was disabled by litc on 2014-07-02;

Usage:
  Options:
  -m <file|dir>  Group of Marker files or directory hold map files,(the filename must be Chr(num)\_\*.map) forced
  -o <file>  Key of output file, forced
  -r <num>	#####行数
  -c <num>	#####列数
  -h         Help

USAGE
	print $usage;
	exit;
}
