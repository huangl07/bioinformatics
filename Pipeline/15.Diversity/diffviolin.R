library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'vcf','v',0,'character',
	'out','o',0,'character',
	'pop','p',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--vcf	the vcf file
	--pop	the group list file 
	--out	the output file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$vcf)){ print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }
if ( is.null(opt$pop)){ print_usage(spec) }

times<-Sys.time()
library(vcfR)
library(ggplot2)
library(reshape2)
vcf<-read.vcfR(file=opt$vcf)
pop<-read.table(opt$pop)
idnames=colnames(vcf@gt)
popid<-NULL
for(i in 2:length(idnames)){popid[i-1]=pop$V2[pop$V1 == idnames[i]]}
vcf<-genetic_diff(vcf,pop=as.factor(popid))
write.table(file=paste(opt$out,".diff.csv",sep=""),row.names=FALSE,vcf)
n=length(unique(pop$V2))
n1=length(colnames(vcf))
n=3+n-1
dpf <- melt(vcf[,c(3:n,n1)], varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
p <- ggplot(dpf, aes(x=variable, y=Depth)) + geom_violin(fill="#2ca25f", adjust = 1.2)
p <- p + xlab("")
p <- p + ylab("")
p <- p + theme_bw()
ggsave(filename=paste(opt$out,".diff.volion.pdf",sep=""),device="pdf",p)
ggsave(filename=paste(opt$out,".diff.volion.png",sep=""),device="png",p)
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
