#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.4.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fSLAF,$fOut,$PNG);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$fSLAF,
				"o:s"=>\$fOut,
				"png"=>\$PNG,
				) or &USAGE;
&USAGE unless ($fSLAF and $fOut);

my $ConvertSVGtoPNG="$Bin/svg2xxx_release/svg2xxx";


# ------------------------------------------------------------------
# Get SLAF info
# ------------------------------------------------------------------
my %SLAF;  my @files;
my $Chro;
if (-f $fSLAF) {
	push @files,$fSLAF;
}elsif(-d $fSLAF){
	@files=glob("$fSLAF/*.map");
}

foreach my $file (@files) {
	open (IN,"<",$file) or die $!;
	($Chro)=$file=~/\/?(\S+)\.(\w+)?\.?map/;

	while (<IN>){
		chomp;
		s/\r//g;
		next if (/^$/ || /^;/ || /^\#/ || /^nloc/) ;
		if (/^group\s+(\w+)/i) {
			next;
		}else{
			my ($Id,$Start)=$_=~/^(\S+)\s+(\S+)/;
			$SLAF{$Chro}{$Id}{"end"}=sprintf("%.3f",$Start);
			$SLAF{$Chro}{$Id}{"id"}=$Id;
			
		}
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


# ------------------------------------------------------------------
# Draw SVG parameter
# ------------------------------------------------------------------
my $LeftMargin=150;
my $RightMargin=10;
my $TopMargin=50;
my $BottomMargin=50;

my $ScaleTxtSpace=45;
my $ScaleIDSpace=60;

my $ScaleLineLen=10;
my $ScaleCrossLine=40;


my $OffsetScale=5;
my $ScaleHeight=15;
my $AnnotationSpace=30;

my $ChroTxtSpace=20;
my $ChroWidth=35;
my $ChroColor="#999900";
my $ChroInterval=35;
my $MiddleColor="rgb(247,251,209)";




# ------------------------------------------------------------------
# Drawing
# ------------------------------------------------------------------
my @SVGFiles=();
foreach my $chro (sort {$a <=> $b} keys %SLAF) {
	my $TotalMarker=scalar keys %{$SLAF{$chro}};



	my $ChroSpace=$ChroWidth+$ScaleTxtSpace+($ScaleLineLen+$ScaleCrossLine)*2+$ScaleIDSpace;
	my $PaperWidth=$LeftMargin+$RightMargin+$ChroSpace+$AnnotationSpace;
	my $PaperHeight=$TopMargin+$ChroTxtSpace+$BottomMargin+$TotalMarker*$ScaleHeight+$ChroWidth;
	my $Scale=$Len{$chro}/($TotalMarker*$ScaleHeight);


	open (SVG,">","$fOut$chro.svg") or die $!;
	push @SVGFiles,"$fOut$chro.svg",$PaperWidth,$PaperHeight;
	print SVG &svg_paper($PaperWidth,$PaperHeight),"\n";

	print SVG "<defs>\n<linearGradient id=\"Gene\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">\n<stop offset=\"0%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n<stop offset=\"50%\" style=\"stop-color:$MiddleColor;stop-opacity:1\"/>\n<stop offset=\"100%\" style=\"stop-color:$ChroColor;stop-opacity:1\"/>\n</linearGradient>\n</defs>","\n";

	my $map_title_prefix=length $chro ==1 ?"LG0":"LG";
	print SVG &svg_txt($LeftMargin,$TopMargin-10,12,"black","$map_title_prefix".$chro),"\n";
	print SVG &svg_circle($LeftMargin+$ChroWidth/2,$TopMargin+$ChroTxtSpace,$ChroWidth/2,$ChroColor),"\n";
	print SVG &svg_circle($LeftMargin+$ChroWidth/2,$TopMargin+$ChroTxtSpace+$Len{$chro}/$Scale,$ChroWidth/2,$ChroColor),"\n";
	print SVG &svg_rect($LeftMargin,$TopMargin+$ChroTxtSpace,$ChroWidth,$Len{$chro}/$Scale,$ChroColor),"\n";

	print SVG &svg_txt($LeftMargin+$ChroWidth+$ScaleLineLen+$ScaleCrossLine+$ScaleIDSpace,$TopMargin+$ChroTxtSpace*2,12,"black","Total Marker:$TotalMarker"),"\n";

	my $numSLAF=1;
	foreach my $start (sort {$SLAF{$chro}{$a}{"end"} <=> $SLAF{$chro}{$b}{"end"}} keys %{$SLAF{$chro}}) {
		my $posi=$LeftMargin;
		print SVG &svg_line($posi-$ScaleLineLen,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{"end"}/$Scale,$posi,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{"end"}/$Scale,"black",1);
		print SVG &svg_line($posi+$ChroWidth,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{"end"}/$Scale,$posi+$ChroWidth+$ScaleLineLen,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{"end"}/$Scale,"black",1);
		print SVG &svg_txt($posi-$ScaleLineLen-$ScaleTxtSpace-$ScaleCrossLine,$TopMargin+$ChroTxtSpace+$numSLAF*$ScaleHeight+$OffsetScale,10,"black",$SLAF{$chro}{$start}{"end"}."cM");
		print SVG &svg_txt($posi+$ChroWidth+$ScaleLineLen+$ScaleCrossLine,$TopMargin+$ChroTxtSpace+$numSLAF*$ScaleHeight+$OffsetScale,10,"black",$SLAF{$chro}{$start}{"id"});
		print SVG &svg_line($posi-$ScaleLineLen,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{'end'}/$Scale,$posi-$ScaleLineLen-$ScaleCrossLine,$TopMargin+$ChroTxtSpace+$numSLAF*$ScaleHeight+$OffsetScale,"black",1);
		print SVG &svg_line($posi+$ChroWidth+$ScaleLineLen,$TopMargin+$ChroTxtSpace+$SLAF{$chro}{$start}{'end'}/$Scale,$posi+$ChroWidth+$ScaleLineLen+$ScaleCrossLine,$TopMargin+$ChroTxtSpace+$numSLAF*$ScaleHeight+$OffsetScale,"black",1);
		$numSLAF++;
	}
	print SVG &svg_end();
}
close (SVG) ;

print STDOUT "Drawing SVG file done.\n";
if ($PNG) {
	for (my $i=0;$i<@SVGFiles ;$i+=3) {
		print STDOUT "Convert $SVGFiles[$i] to PNG format ... \n";
		my $dirname=dirname($SVGFiles[$i]);
		my $basename=basename($SVGFiles[$i]);
		`mkdir $dirname` unless (-d $dirname) ;
		my $pwd=`pwd`;chomp $pwd;
		chdir $dirname;
		`$ConvertSVGtoPNG $basename -w $SVGFiles[$i+1] -h $SVGFiles[$i+2]`;
		chdir $pwd;

	}
}


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
Program: SLAFonRef.pl (遗传图谱绘图程序)
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>
Description：
	Support JoinMap‘s map format files;
	Support Converting svg to png files;
	Every linkage group will be drawn in each file;
	Auto OffSet was disabled by litc on 2011-04-26;

Usage:
  Options:
  -m <file|dir>  Group of Marker files or directory hold map files, forced
  -o <file>  Key of output file, forced
  -png       Convert svg to png, default off
  -h         Help

USAGE
	print $usage;
	exit;
}
