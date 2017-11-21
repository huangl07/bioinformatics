#!/usr/bin/env Rscript

times<-Sys.time()
library('getopt');
options(bitmapType='cairo')

spec = matrix(c(
	'i','a',0,'character',
	'o','b',0,'character',
	'help','c',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript insert_size.r --i  --o  
	
Usage:
	--i     insert_size file
	--o	    insert_size picture name
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }
a=read.table(opt$i)
pdffile=paste(opt$o,".pdf",sep="");
pdf(file=pdffile)
len<-length(a$V2);
if (len < 1000){
	ymax<-max(a$V2[0:len]);
}else{
	ymax<-max(a$V2[0:1000]);
}
plot(a[,1],a[,2],type='h',col="#0099FF",xlab="Insert Size(bp)",ylab="Reads Number",main="Insert Size Distribution",xlim=c(0,700),ylim=c(0,ymax*1.2))
dev.off()
pdffile=paste(opt$o,".png",sep="");
png(file=pdffile)
plot(a[,1],a[,2],type='h',col="#0099FF",xlab="Insert Size(bp)",ylab="Reads Number",main="Insert Size Distribution",xlim=c(0,700),ylim=c(0,ymax*1.2))
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
#TSHKO.insertsize.pdf
