library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	the input  file
	--output	the out file 
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

d<-read.table(opt$input,head=TRUE,comment.char="^")
d$pvalue=1-d$pvalue;
d$Nmiss=1-d$Nmiss/d$Nind;

library(qqman)
pdf(paste(opt$out,".seg.pdf"),height=900,width=1600)
manhattan(d,chr="GroupID",bp="Dis",p="pvalue",logp=FALSE,col=rainbow(4),suggestiveline=0.95,genomewideline=FALSE,ylab="Segregation Pvalue",main="Marker Disortion Distribution")
box()
dev.off()
png(paste(opt$out,".seg.png"),height=900,width=1600)
manhattan(d,chr="GroupID",bp="Dis",p="pvalue",logp=FALSE,col=rainbow(4),suggestiveline=0.95,genomewideline=FALSE,ylab="Segregation Pvalue",main="Marker Disortion Distribution")
box()
dev.off()
chr=unique(d$GroupID);
for(i in chr){
	df=data.frame(chr=d$GroupID[d$GroupID==i],bp=d$Dis[d$GroupID==i],p=d$pvalue[d$GroupID==i]);
	pdf(paste(opt$out,".LG",i,".seg.pdf",sep=""),height=900,width=1600)
	manhattan(df,chr="chr",bp="bp",p="p",logp=FALSE,col="blue",suggestiveline=FALSE,genomewideline=0.95,ylab="Segregation Pvalue",type="l",main="Marker Disortion Distribution")
	box()
	dev.off()
	png(paste(opt$out,".LG",i,".seg.png",sep=""),height=900,width=1600)
	manhattan(df,chr="chr",bp="bp",p="p",logp=FALSE,col="blue",suggestiveline=FALSE,genomewideline=0.95,ylab="Segregation Pvalue",type="l",main="Marker Disortion Distribution")
	box()
	dev.off()

}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
