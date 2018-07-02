library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'thre','t',1,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	 the input  file
	--output	the out file 
	--thre    threshold
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

bsa<-read.table(opt$input,head=TRUE)
library(qqman)
thre<-as.integer(opt$thre)
pdf(paste(opt$output,"pdf",sep="."),height=9,width=16)
manhattan(bsa,chr="CHROM",bp="POS",p="Gprime",col=rainbow(4),logp=FALSE,ylab="chrosome",xlab="G value",suggestiveline=thre,cex=2)
dev.off()
png(paste(opt$output,"png",sep="."),height=900,width=1600)
manhattan(bsa,chr="CHROM",bp="POS",p="Gprime",col=rainbow(4),logp=FALSE,ylab="Genomic Position()",xlab="G value",suggestiveline=thre,cex=2)
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)