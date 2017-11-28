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

GWAS<-read.table(opt$input,head=TRUE)
len=length(GWAS$p_value);
GWAS$p_value[GWAS$p_value == 0]= 0.0001/len;
print(0.0001/len)
library(qqman)
pdf(paste(opt$output,".manhattan.pdf",sep=""),height=900,width=1600)
manhattan(GWAS,chr="chr",bp="pos",p="p_value",logp=TRUE,col=rainbow(4),suggestiveline=FALSE,genomewideline=0.0001/len)
dev.off()
png(paste(opt$output,".manhattan.png",sep=""),height=900,width=1600)
manhattan(GWAS,chr="chr",bp="pos",p="p_value",logp=TRUE,col=rainbow(4),suggestiveline=FALSE,genomewideline=0.0001/len)
dev.off()
pdf(paste(opt$output,".qq-plot.pdf",sep=""))
qq(GWAS$p_value)
dev.off()
png(paste(opt$output,".qq-plot.png",sep=""))
qq(GWAS$p_value)
dev.off()

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
