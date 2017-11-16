library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'sliding','s',0,'character',
	'output','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--sliding	the input sliding file
	--output	the out file 
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$sliding)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

sliding<-read.table(opt$sliding,stringsAsFactors=FALSE)
library(qqman)
names(sliding)<-c("chr","pos1","pos2","index1","index2","index3","snpnum")
chrlab=unique(sliding$chr)
for (i in (1:length(chrlab))){sliding$chr[which(sliding$chr==chrlab[i])]=i}
sliding$chr=as.numeric(sliding$chr);
pdf(paste(opt$output,".pdf",sep=""),height=900,width=1600)
par(mfrow = c(3, 1))
manhattan(sliding,chr="chr",bp="pos2",p="index1",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index 1",ylim=c(0,1),logp=FALSE)
manhattan(sliding,chr="chr",bp="pos2",p="index2",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index 2",ylim=c(0,1),logp=FALSE)
manhattan(sliding,chr="chr",bp="pos2",p="index3",col=rainbow(4),chrlabs=chrlab,ylab="delta SNP-index",ylim=c(-1,1),logp=FALSE)
dev.off()
png(paste(opt$output,".png",sep=""),height=900,width=1600)
par(mfrow = c(3, 1))
manhattan(sliding,chr="chr",bp="pos2",p="index1",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index 1",ylim=c(0,1),logp=FALSE)
manhattan(sliding,chr="chr",bp="pos2",p="index2",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index 2",ylim=c(0,1),logp=FALSE)
manhattan(sliding,chr="chr",bp="pos2",p="index3",col=rainbow(4),chrlabs=chrlab,ylab="delta SNP-index",ylim=c(-1,1),logp=FALSE)
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
