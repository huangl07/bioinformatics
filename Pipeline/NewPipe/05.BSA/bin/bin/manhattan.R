library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'col','c',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	the input file
	--output	the out file \
	--col	the col number for draw
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

data<-read.table(opt$input,stringsAsFactors=FALSE,head=FALSE)
data<-data[-1,]
library(qqman)
collist <- unlist(strsplit(opt$col, split=","))
if(length(collist) ==4){
	chr<-data[[as.numeric(collist[1])]]
	pos<-data[[as.numeric(collist[2])]]
	index<-data[[as.numeric(collist[3])]]
	thres<-data[[as.numeric(collist[4])]]
	chrlab=unique(chr);
	for (i in 1:length(chrlab)){chr[chr==chrlab[i]]=i}
	df<-data.frame(chr=as.numeric(chr),pos=as.numeric(pos),index=as.numeric(index))
	pdf(paste(opt$output,".pdf",sep=""),height=900,width=1600)
	manhattan(df,chr="chr",bp="pos",p="index",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index",ylim=c(0,1),logp=FALSE,suggestiveline=unique(as.numeric(thres)),main="SNP INDEX DISTRIBUTION",type="h",)
	dev.off()
	png(paste(opt$output,".png",sep=""),height=900,width=1600)
	manhattan(df,chr="chr",bp="pos",p="index",col=rainbow(4),chrlabs=chrlab,ylab="SNP-index",ylim=c(0,1),logp=FALSE,suggestiveline=unique(as.numeric(thres)),main="SNP INDEX DISTRIBUTION",type="h")
	dev.off()
}else{
	chr<-data[[as.numeric(collist[1])]]
	pos<-data[[as.numeric(collist[2])]]
	index1<-data[[as.numeric(collist[3])]]
	index2<-data[[as.numeric(collist[4])]]
	delta<-data[[as.numeric(collist[5])]]
	chrlab=unique(chr);
	for (i in 1:length(chrlab)){chr[chr==chrlab[i]]=i}
	df<-data.frame(chr=as.numeric(chr),pos=as.numeric(pos),index1=as.numeric(index1),index2=as.numeric(index2),delta=as.numeric(delta))
	pdf(paste(opt$output,".pdf",sep=""),height=900,width=1600)
	par(mfrow = c(3, 1))
	manhattan(df,chr="chr",bp="pos",p="index1",col=rainbow(4),chrlabs=chrlab,ylab="Bulk SNP-index 1",ylim=c(0,1),logp=FALSE)
	manhattan(df,chr="chr",bp="pos",p="index2",col=rainbow(4),chrlabs=chrlab,ylab="Bulk SNP-index 2",ylim=c(0,1),logp=FALSE)
	manhattan(df,chr="chr",bp="pos",p="delta",col=rainbow(4),chrlabs=chrlab,ylab="Delta SNP-index",ylim=c(-1,1),logp=FALSE,suggestiveline=unique(as.numeric(thres)))
	dev.off()
	png(paste(opt$output,".png",sep=""),height=900,width=1600)
	par(mfrow = c(3, 1))
	manhattan(df,chr="chr",bp="pos",p="index1",col=rainbow(4),chrlabs=chrlab,ylab="Bulk SNP-index 1",ylim=c(0,1),logp=FALSE)
	manhattan(df,chr="chr",bp="pos",p="index2",col=rainbow(4),chrlabs=chrlab,ylab="Bulk SNP-index 2",ylim=c(0,1),logp=FALSE)
	manhattan(df,chr="chr",bp="pos",p="delta",col=rainbow(4),chrlabs=chrlab,ylab="Delta SNP-index",ylim=c(-1,1),logp=FALSE,suggestiveline=unique(as.numeric(thres)))
	dev.off()
}
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
