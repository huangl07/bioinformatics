times<-Sys.time() 
library(qqman)
library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
    'infile','f',0,'character',
    'outfile','o',0,'character',
    'table','l',0,'character'
     ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
    Usage example: 

    Usage:
    --infile sweep out file
    --outfile figure name
    --table out file  
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }

mydata<-read.table(opt$infile,header = TRUE,sep="\t")
ymaxP=ceiling(max(mydata$P,na.rm=TRUE))
png(paste(opt$outfile,".png",sep=""),height=900,width=1600)
manhattan(mydata,chr="CHR",bp="POS",p="P",col=rainbow(4),chrlabs=("SNPs"),ylab="log10(BF)",ylim=c(0,ymaxP),suggestiveline=FALSE ,genomewideline=FALSE,logp=FALSE)
dev.off()
pdf(paste(opt$outfile,".pdf",sep=""),height=8,width=16)
manhattan(mydata,chr="CHR",bp="POS",p="P",col=rainbow(4),chrlabs=("SNPs"),ylab="log10(BF)",ylim=c(0,ymaxP),suggestiveline=FALSE,genomewideline=FALSE,logp=FALSE)
dev.off()
select <- subset(mydata,log10(mydata$P)>=1.5)
select <- select[with(select,order(mydata$POS)),]
select <- select[order(mydata$POS,mydata$P),]
write.table(select, file = opt$table)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
