#! /mnt/ilustre/app/pub/R/bin/Rscript

#Args <- commandArgs()
#m<-read.table(Args[6],header=T)
#b<-read.table(Args[7],header=T)

library(getopt)
options(bitmapType='cairo')

opt = getopt(matrix(c(
  'input','i',1,'character',
  'output','o',2,'character',
  'help','h',0,'logical'
),byrow=TRUE, ncol=4));
usage<-function(){
cat("This script is used to draw colinear figure
Usage	Rscript colinearity-plot.R for fig 3-14 [options]
Options:
	-input, --input	mapinfo file
	-output, --output file name
	-h, --help	print display this help and exit
")
  q(status=1);
}
if (!is.null(opt$help) ) { usage() }
if (is.null(opt$input) ) { usage() }
if (is.null(opt$output) ) { usage() }

mb<-read.table(opt$input,header=T)
library(MASS)
library(plyr)
mbs<-mb[order(mb$chr),]
x<-ddply(mbs,"chr",summarise,max=max(pos),min=min(pos))
num<-c(0,x$max[1:length(x$max)-1])-x$min
d1<-data.frame()
for(i in 1:length(x$chr)){
s=data.frame()
s=subset(mbs,chr==x$chr[i])
s$pos=s$pos+sum(num[1:i])
d1=rbind(d1,s)
}
d2<-d1[order(d1$pos),]

st<-ddply(d2,c("chr","lg"),summarise,num=table(lg))
k<-ddply(st,"lg",summarise,max=max(num))
st$int<-interaction(st$lg,st$num)
k$int<-interaction(k$lg,k$max)
re<-st[match(k$int,st$int),]

y<-ddply(d2,"lg",summarise,max=max(cm),min=min(cm),meam=mean(pos))
rf<-merge(re[order(re$chr),],y,by="lg")
rg<-rf[order(rf$meam),]
nu<-c(0,rg$max[1:length(rg$max)-1])-rg$min
d3<-data.frame()
for(i in 1:length(rg$lg)){
s=data.frame()
s=subset(d2,lg==rg$lg[i])
s$cm=s$cm+sum(nu[1:i])
d3=rbind(d3,s)
}
d4<-d3[order(d3$pos),]

chr<-ddply(d4,"chr",summarise,pos=(max(pos)+min(pos))/2)
lg<-ddply(d4,c("chr","lg"),summarise,cm=(max(cm)+min(cm))/2)
lg$int<-interaction(lg$chr,lg$lg)
re$int<-interaction(re$chr,re$lg)
lg2<-lg[match(re$int,lg$int),]
lg3<-lg2[order(lg2$chr),]
v<-ddply(d4,"chr",summarise,pos=max(pos))
h<-ddply(d4,c("chr","lg"),summarise,cm=max(cm))
h$int<-interaction(h$chr,h$lg)
h1<-h[match(re$int,h$int),]
h2<-h1[order(h1$chr),]

library(ggplot2)
chrs<-as.character(d4$chr)
col<-colors()[c(rep(c(26,258,24,552),times=ceiling(length(unique(chrs))/4))[1:length(unique(chrs))])]
p=ggplot(d4,aes(x=pos,y=cm,colour=chrs))+geom_point(shape=20,size=0.5)+scale_colour_manual(values=col,guide=FALSE)+xlab("Chromosome") + ylab("Linkage Group")+scale_y_continuous(breaks=lg3$cm,labels=lg3$lg)+scale_x_continuous(breaks=chr$pos,labels=chr$chr)+theme_bw()+labs(title="Collinary by Reference")
p=p+theme(panel.grid.major = element_blank(),plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank())+geom_vline(xintercept=v$pos[1:length(v$pos)-1],size=0.1,colour="gray80")+geom_hline(yintercept=h2$cm[1:length(h2$cm)-1],size=0.1,colour="gray80")
lgfac<-factor(d4$lg,ordered=T)
table318<-data.frame()
for (j in levels(lgfac)){
	d4tmp=d4[which(d4$lg==j),]
	spcor<-cor(d4tmp$pos,d4tmp$cm,method ="spearman")
	table318[j,1]=j
	table318[j,2]=spcor
}
colnames(table318)=c("LG ID","Spearman")

ggsave(paste(opt$output,".pdf",sep=""), p)
ggsave(paste(opt$output,".png",sep=""), p)

write.table(table318,file=paste(opt$output,".xls",sep=""),row.names=F,sep="\t")
