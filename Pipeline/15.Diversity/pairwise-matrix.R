library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'genepop','g',0,'character',
	'output','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--genepop	the input file
	--output	the out file \
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$genepop)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()
library("FinePop")
popdata <- read.genepop(genepop=opt$genepop)
draw.matrix<-function(res=NULL,out=NULL,type=NULL){
	poplevels=levels(as.factor(rownames(res)))
	numericMatrix=res
	Matrix <- numericMatrix
	ColorRamp <- colorRampPalette(c("white", "steelblue1", "blue3"))
	outfileGraphic <- paste(out,type,"png", sep=".")
	png(outfileGraphic, width=1300, height=1300, res=144)
	smallplot <- c(0.874, 0.9, 0.18, 0.83)
	bigplot <- c(0.13, 0.85, 0.14, 0.87)
	old.par <- par(no.readonly = TRUE)
	par(plt = smallplot)
	Min <- min(Matrix, na.rm=TRUE)
	Max <- max(Matrix, na.rm=TRUE)
	binwidth <- (Max - Min) / 64
	y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
	z <- matrix(y, nrow = 1, ncol = length(y))
	image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)
	if(Min == Max){
		axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
	} else {
		axis(side=4, las = 2, cex.axis=0.8)
	}
	box()
	mtext(text=expression(bold(Diversity~Parameter)), side=4, line=2.5, cex=1.1)
	a <- ncol(numericMatrix)
	b <- nrow(numericMatrix)
	x <- c(1:a)
	y <- c(1:b)
	par(new = TRUE, plt = bigplot)
	image(x,y,as.matrix(Matrix), col=ColorRamp(64),
			main=expression(bold(Matrix~of~pairwise~diversity)), xlab="",
			ylab="", axes=FALSE)
	box()
	Labels=poplevels
	if(is.null(Labels)){
		axis(1, at = c(1:a))
		axis(2, at = c(1:b), labels=c(b:1))
		mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
			font=2)
		mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
			font=2)
	} else{
		axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
			 las=2)
		axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
			 las=2)
	}
	par(old.par)#reset graphic parameters
	dev.off()
	outfileGraphic <- paste(out,type,"pdf", sep=".")
	pdf(outfileGraphic, width=13, height=13)
	smallplot <- c(0.874, 0.9, 0.18, 0.83)
	bigplot <- c(0.13, 0.85, 0.14, 0.87)
	old.par <- par(no.readonly = TRUE)
	par(plt = smallplot)
	Min <- min(Matrix, na.rm=TRUE)
	Max <- max(Matrix, na.rm=TRUE)
	binwidth <- (Max - Min) / 64
	y <- seq(Min + binwidth/2, Max - binwidth/2, by = binwidth)
	z <- matrix(y, nrow = 1, ncol = length(y))
	image(1, y, z, col = ColorRamp(64),xlab="", ylab="", axes=FALSE)
	if(Min == Max){
		axis(side=4, las = 2, cex.axis=0.8, at=Min, labels=round(Min, 2))
	} else {
		axis(side=4, las = 2, cex.axis=0.8)
	}
	box()
	mtext(text=expression(bold(Diversity~Parameter)), side=4, line=2.5, cex=1.1)
	a <- ncol(numericMatrix)
	b <- nrow(numericMatrix)
	x <- c(1:a)
	y <- c(1:b)
	par(new = TRUE, plt = bigplot)
	image(x,y,as.matrix(Matrix), col=ColorRamp(64),
			main=expression(bold(Matrix~of~pairwise~diversity)), xlab="",
			ylab="", axes=FALSE)
	box()
	Labels=poplevels
	if(is.null(Labels)){
		axis(1, at = c(1:a))
		axis(2, at = c(1:b), labels=c(b:1))
		mtext(side = 1, at =(a/2), line = 2.5, text = "Population", cex=1,
			font=2)
		mtext(side = 2, at =(b/2), line = 2.7, text = "Population", cex=1,
			font=2)
	} else{
		axis(1, at = c(1:a), labels=Labels[1:length(Labels)], cex.axis=0.75,
			 las=2)
		axis(2, at = c(1:b), labels=Labels[length(Labels):1], cex.axis=0.75,
			 las=2)
	}
	par(old.par)#reset graphic parameters
	dev.off()
}

GstH=GstH(popdata)
GstN=GstN(popdata);
GstNC=GstNC(popdata)
Fst=thetaWC.pair(popdata)
DJ=DJ(popdata)
write.table(file=paste(opt$output,"GstN.matrix",sep="."),GstN,sep="\t")
write.table(file=paste(opt$output,"Fst.matrix",sep="."),Fst,sep="\t")
write.table(file=paste(opt$output,"GstNC.matrix",sep="."),GstNC,sep="\t")
write.table(file=paste(opt$output,"DJ.matrix",sep="."),DJ,sep="\t")
write.table(file=paste(opt$output,"GstH.matrix",sep="."),GstH,sep="\t")
draw.matrix(res=GstN,out=opt$output,type="GstN")
draw.matrix(res=Fst,out=opt$output,type="Fst")
draw.matrix(res=GstNC,out=opt$output,type="GstNC")
draw.matrix(res=DJ,out=opt$output,type="DJ")
draw.matrix(res=GstH,out=opt$output,type="GstH")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
