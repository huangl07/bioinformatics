library(qqman)
options(scipen=200)
snp<-read.table("index-calc.result.index",comment.char="^",stringsAsFactors=F,head=T)
win<-read.table("sliding-win.result",comment.char="^",stringsAsFactors=F,head=T)
chr<-snp$X.chr
pos<-snp$pos
index<-snp$INDEX2
chrlab=unique(chr);
for (i in 1:length(chrlab)){chr[chr==chrlab[i]]=i}
dfsnp<-data.frame(chr=as.numeric(chr),pos=as.numeric(pos),index=as.numeric(index))

chr<-win$chr
pos<-win$pos1
index<-win$index2
chrlab=unique(chr);
for (i in 1:length(chrlab)){chr[chr==chrlab[i]]=i}
dfwin<-data.frame(chr=as.numeric(chr),pos=as.numeric(pos),index=as.numeric(index))
for (i in 1:length(chrlab)){
	pdf(paste(chrlab[i],"pdf",sep="."))
	newdfsnp<-subset(dfsnp,dfsnp$chr == i);
	newdfwin<-subset(dfwin,dfwin$chr == i);
	manhattan(newdfsnp,chr="chr",bp="pos",p="index",col="green",type="p",ylab="SNP-index",xlab="Chromesome Wide Distribution(kb)",ylim=c(0,1),logp=FALSE,main=paste(chrlab[i],"SNP INDEX DISTRIBUTION",sep=" "))
	lines(newdfwin$pos,newdfwin$index,type="l",lwd=2)
	dev.off()
	png(paste(chrlab[i],"png",sep="."))
	manhattan(newdfsnp,chr="chr",bp="pos",p="index",col="green",type="p",ylab="SNP-index",xlab="Chromesome Wide Distribution(kb)",ylim=c(0,1),logp=FALSE,main=paste(chrlab[i],"SNP INDEX DISTRIBUTION",sep=" "))
	lines(newdfwin$pos,newdfwin$index,type="l",lwd=2)
	dev.off()

}