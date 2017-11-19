library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'infile','f',0,'character',
	'outfile','o',0,'character',
	'type','t',0,'character',
	'lab','l',0,'character',
	'str','s',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--infile	the input omega file
	--type	the input sweed file
	--outfile	the out file 
	--lab	the ylab
	--str	the label for draw eg(chr,pos,omega)
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help)) { print_usage(spec) }
if ( is.null(opt$infile)) { print_usage(spec)}
if ( is.null(opt$type)){ opt$type="p" }
if ( is.null(opt$str)){ print_usage(spec) }

times<-Sys.time()
infile<-read.table(opt$infile,stringsAsFactors=FALSE,head=TRUE,na.strings="nan")
label=strsplit(opt$str,split=",")
if ( is.null(opt$lab)){ opt$lab=label[[1]][3]}

library(qqman)
chr=infile[[label[[1]][1]]]
pos=infile[[label[[1]][2]]]
pvalue=infile[[label[[1]][3]]]
draw<-data.frame(CHR=chr,POS=pos,P=pvalue)

chrlab1=unique(draw$CHR)
for (i in (1:length(chrlab1))){draw$CHR[which(draw$chr==chrlab1[i])]=i}
ymaxP=ceiling(max(draw$P)/0.8)
pdf(paste(opt$outfile,".pdf",sep=""),height=900,width=1600)
manhattan(draw,chr="CHR",bp="POS",p="P",col=rainbow(4),chrlabs=chrlab1,ylab=opt$lab,ylim=c(0,ymaxP),type=opt$type,suggestiveline=FALSE,genomewideline=FALSE,logp=FALSE)
dev.off()
png(paste(opt$outfile,".png",sep=""),height=900,width=1600)
manhattan(draw,chr="CHR",bp="POS",p="P",col=rainbow(4),chrlabs=chrlab1,ylab=opt$lab,ylim=c(0,ymaxP),type=opt$type,suggestiveline=FALSE,genomewideline=FALSE,logp=FALSE)
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
