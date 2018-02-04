library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','f',0,'character',
	'out','o',0,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--input	the raisd file
	--out the output file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }
if ( is.null(opt$input)){ print_usage(spec) }

library(qqman)
times<-Sys.time()
data<-read.table(opt$input,sep="\t",head=TRUE,stringsAsFactors=FALSE)
thres<-quantile(data$MuStat,probs=0.999)
write.table(file=paste(opt$out,".raisd",sep=""),row.names=FALSE,data)
write.table(file=paste(opt$out,".raisd.select",sep=""),row.names=FALSE,subset.data.frame(data,data$MuStat >= thres))
chrlab=unique(data$chr)
for (i in 1:length(chrlab)){data$chr[data$chr==chrlab[i]]=i}
data$chr<-as.numeric(data$chr)

pdf(paste(opt$out,".raisd.pdf",sep=""))
manhattan(data,chr="chr",bp="pos",p="MuStat",col=rainbow(4),chrlabs=chrlab,ylab="Pi",logp=FALSE,suggestiveline=FALSE,ylim=c(0,max(data$MuStat)/0.8))
box()
dev.off()
png(paste(opt$out,".raisd.png",sep=""))
manhattan(data,chr="chr",bp="pos",p="MuStat",col=rainbow(4),chrlabs=chrlab,ylab="Pi",logp=FALSE,suggestiveline=FALSE,ylim=c(0,max(data$MuStat)/0.8))
box()
dev.off()


escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
