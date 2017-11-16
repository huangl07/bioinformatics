#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'help','h',0,"logical",
	'infile','i',1,"character",
	'outfile','o',1,"character"
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
Options: 
	--help		NULL 		get this help
	--infile 	character 	structure file
	--outfile 	character 	the filename for output graph [forced]
	")
	q(status=1);
}
if(is.null(opt$infile)){print_usage(spec)}
if(is.null(opt$outfile)){print_usage(spec)}

tbl<-read.table(opt$infile,header=FALSE);
pdf(paste(opt$outfile,".pdf",sep=""))
plot(tbl$V1,tbl$V2,xlab="K value",ylab="CV error",type="b",lty=1,col="blue")
axis(tbl$V1,at=tbl$V1)
points(x=tbl$V1[which.min(tbl$V2)],y=min(tbl$V2),pch=25, col="red", bg="blue")
dev.off()
png(paste(opt$outfile,".png",sep=""))
plot(tbl$V1,tbl$V2,xlab="K value",ylab="CV error",type="b",lty=1,col="blue")
axis(tbl$V1,at=tbl$V1)
points(x=tbl$V1[which.min(tbl$V2)],y=min(tbl$V2),col="red")
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
