#!/usr/bin/env Rscript

times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',0,'character',
	'output','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript indel_len.r --i  --o  
	
Usage:
	--i	insert_size file
	--o	the output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input) ) { print_usage(spec) }
if ( is.null(opt$output) ) { print_usage(spec) }
df<-read.table(opt$input,header=FALSE)
col=length(colnames(df))
suml<-function(x,df){
	a<-df[df$V1 ==x,2:col]
	return(sum(a))
}
x=data.frame(id=df$V1);
twin=apply(x,MARGIN=1,function(x,y) suml(x[1],df));
tbl=data.frame(id=df$V1,num=twin);
print(tbl);
pdf(paste(opt$output,".pdf",sep=""))
plot(tbl,type="l",main="length distribution",col="blue",xlim=c(0,5000))
dev.off()
png(paste(opt$output,".png",sep=""))
plot(tbl,type="l",main="length distribution",col="blue",xlim=c(0,5000))
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)





#TSHKO.insertsize.pdf
