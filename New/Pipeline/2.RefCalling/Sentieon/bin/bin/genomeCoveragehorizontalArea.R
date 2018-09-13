#!/usr/bin/env Rscript


# load library
library('getopt');
library(ggplot2);
library(methods);
options(bitmapType='cairo');


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'help' , 'h', 0, "logical",
	'infile' , 'i', 1, "character",
	'idfile' , 'd', 1, "character",
	'outfile' , 'o', 1, "character",
	'group.col' , 'g', 1, "integer",
	'x.col' , 'x', 1, "integer",
	'y.col' , 'y', 1, "integer",
	'height' , 'H', 1, "integer",
	'width' , 'W', 1, "integer",
	'x.lab' , 'X', 1, "character",
	'y.lab' , 'Y', 1, "character",
	'title.lab' , 'T', 1, "character",
	'unit' , 'u', 1, "character",
	'lab.size' , 'l', 1, "integer",
	'strip.size' , 't', 1, "integer",
	'axis.x.size' , 's', 1, "integer",
	'axis.y.size' , 'z', 1, "integer",
	'no.grid' , 'r', 0, "logical",
	'log2','2',0,"logical",
	'skip' , 'k', 1, "integer"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
1) Rscript genomeCoveragehorizontalArea.r --infile in_genome_cov.data --idfile chr.txt \\
	--outfile out_chr_cov.png --group.col 1 --x.col 2 --y.col 3 \\
	--x.lab \"Chromsome\" --y.lab \"y lab\" --skip 1 --unit \"Mb\"
2) Rscript genomeCoveragehorizontalArea.r --infile in_genome_cov.data --idfile chr.txt \\
	--outfile out_chr_cov.png --group.col 1 --x.col 2 --y.col 3 \\
	--x.lab \"Chromsome\" --y.lab \"y lab\" --skip 1 --axis.y.size 3 --strip.size 12

Options: 
--help		-h 	NULL 		get this help
--infile 	-i 	character 	the input file [forced]
--idfile 	-d 	character 	the input id file [forced]
--outfile 	-o 	character 	the filename for output graph [forced]
--group.col 	-g 	integer 	the col for group factor [forced]
--x.col 	-x 	integer 	the col for x value [forced]
--y.col 	-y 	integer 	the col for y value [forced]
--height 	-H 	integer 	the height of graph [optional, default: 3000]
--width 	-W 	integer 	the width of graph [optional, default: 4000]
--x.lab 	-X 	character 	the lab for x [forced]
--y.lab 	-Y 	character 	the lab for y [forced]
--title.lab 	-T 	character 	the lab for title [optional, default: NULL]
--unit 		-u 	character 	the lab unit for x axis [optional, default: kb]
--lab.size 	-l 	integer 	the font size of lab [optional, default: 10]
--strip.size 	-t 	integer 	the font size of strip lab [optional, default: 8]
--axis.x.size 	-s 	integer 	the font size of text for axis x [optional, default: 8]
--axis.y.size 	-z 	integer 	the font size of text for axis y [optional, default: 5]
--no.grid	-r 	NULL 		Do not drawing grid
--skip 		-k 	integer 	the number of line for skipping [optional, default: 0]
--log2	-2	NULL	log2 or not
\n")
	q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$idfile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ opt$x.col=2 }
if ( is.null(opt$y.col) )	{ opt$y.col=3 }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 10 }
if ( is.null(opt$strip.size ) )		{ opt$strip.size = 8 }
if ( is.null(opt$axis.x.size ) )	{ opt$axis.x.size = 8 }
if ( is.null(opt$axis.y.size ) )	{ opt$axis.y.size = 5 }
if ( is.null(opt$title.lab) )		{ opt$title.lab = NULL }
if ( is.null(opt$unit) )		{ opt$unit = "kb" }
if ( is.null(opt$group.col) )	{ opt$group.col=1 }



#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip,header=F)
# check dim
if (!is.null(opt$log2)){data$V3=log2(data$V3)}

data.dim <- dim(data)
if ( is.null(data.dim) ){
	cat("Final Error: the format of infile is error, dim(data) is NULL\n")
	print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col, opt$group.col) ){
	cat("Final Error: max(x.col, y.col, group.col) > the col of infile\n")
	print_usage(spec)
}
# reading id
id <- read.table(opt$idfile)
id <- as.matrix(id[,1])
# check dim
id.dim <- dim(id)
if ( is.null(id.dim) ){
	cat("Final Error: the format of idfile is error, dim(data) is NULL\n")
	q(status=1)
}
# check col size
if ( id.dim[2] != 1 ){
	cat("Final Error: the col of idfile != 1\n")
	q(status=1)
}
# check and extract
if( length(id[,1]) <= 1 ){
	#cat("Final Error: the row of idfile <= 1\n")
	#q(status=1)
}
if( sum(id[,1] %in% data[,opt$group.col]) != length(id[,1]) ){
	cat("Final Error: id not exists in the infile\n")
	q(status=1)
}
# create df
data <- data[data[,opt$group.col] %in% id[,1],]
df <- data.frame(x=data[,opt$x.col]-1, group=factor(data[,opt$group.col], levels=id[,1]), y=data[,opt$y.col])



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
p <- ggplot(df, aes(x=x, y=y),binwidth=0.1) + geom_area(aes(color=group, fill=group))
p <- p + facet_grid(group~.)
p <- p + xlim(min(df$x), max(df$x))


#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
# lab
p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) + labs(title=opt$title.lab)
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
	axis.text.x = element_text(face="bold", size=opt$axis.x.size),
	axis.text.y  = element_text(face="bold", size=opt$axis.y.size),
	strip.text.y = element_text(size=opt$strip.size, angle=0),
	strip.background = element_rect(colour="white", fill="white"),
	axis.line.x=element_line(),
        axis.line.y=element_line())
# x axis label
x.lab.fun <- function(x){
	lab <- paste(x,"(", opt$unit,")", sep="")
}
p <- p + scale_x_continuous(label=x.lab.fun)
# remove legend
p <- p + theme(legend.position = "none")
# grid and background
if ( !is.null(opt$no.grid) ) {
	#p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
	p <- p + theme( panel.background = element_rect(colour="white", size=1, fill="white"),
		panel.grid = element_blank())
}

#add theme
p<-p+theme(strip.background=element_rect(fill="white", colour="white"), strip.text.y=element_text(angle=360, vjust=0.1), 
		panel.background=element_rect(fill="white", size=10), panel.grid=element_blank(), 
		axis.line=element_line(colour="black"),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))
#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=paste(opt$outfile,".png",sep=""), height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







