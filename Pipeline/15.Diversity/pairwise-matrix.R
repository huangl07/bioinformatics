library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'genepop','g',0,'character',
	'output','o',0,'character',
	'poplist','p',0,'character',
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
	--poplist	the col number for draw
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$genepop)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$poplist)){ print_usage(spec) }
times<-Sys.time()
library("FinePop")
popdata <- read.genepop(genepop=opt$genepop)
GstH=GstH(popdata)
GstN=GstN(popdata);
GstNC=GstNC(popdata)
Fst=thetaWC.pair(popdata)
DJ=DJ(popdata)
write.csv(file=paste(opt$output,"GstN.matrix",sep="."),GstN,sep="\t")
write.csv(file=paste(opt$output,"Fst.matrix",sep="."),GstN,sep="\t")
write.csv(file=paste(opt$output,"GstNC.matrix",sep="."),GstNC,sep="\t")
write.csv(file=paste(opt$output,"DJ.matrix",sep="."),DJ,sep="\t")
write.csv(file=paste(opt$output,"GstH.matrix",sep="."),GstH,sep="\t")

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
