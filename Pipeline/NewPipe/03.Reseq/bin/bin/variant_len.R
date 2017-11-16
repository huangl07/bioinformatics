times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'i','a',0,'character',
	'o','b',0,'character',
	'c','c',0,'character',
	'help','c',0,'logical'
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
	--c	the col string
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$i) ) { print_usage(spec) }
if ( is.null(opt$o) ) { print_usage(spec) }
df<-read.table(opt$i,header=FALSE)
id<-df$V1;
col=length(colnames(df))
tbl<-df[,2:col]/sum(df[,2:col])
pdf(paste(opt$o,"pdf",sep=""))
barplot(t(as.matrix(tbl)),names.arg=id, col=rainbow(col),xlab="length", ylab="percentage", border=NA,las=2)
dev.off()
png(paste(opt$o,"png",sep=""))
barplot(t(as.matrix(tbl)),names.arg=id, col=rainbow(col),xlab="length", ylab="percentage", border=NA,las=2)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)





#TSHKO.insertsize.pdf
