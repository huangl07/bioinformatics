#!/usr/bin/env Rscript
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'hb','p',1,'character',
	'lb','l',1,'character',
	'pop','t',2,'character',
	'bs','b',2,'character',
	'ws','w',2,'character',
	'pta','c',2,'character',
	'ptb','d',2,'character',
	'output','o',1,'character',
	'opt','u','2','character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example:
	Rscript G.R --input --hb --lb --output --all --opt
Usage:
	--input *vcf.table (the vcf to TABLE from GATK table)
	--chr the number of chromosomes
	--hb	HighBulk sampleID
	--lb	LowBulk sampleID
	--pop   all kinds of population
	--bs    bulk size
	--output	output dir
	--ws window size
	--opt save the file name of the p-value less 0.01
	--help		usage
\n")
	q(status=1);
}

times<-Sys.time()

if (is.null(opt$input)) { print_usage(spec)}
if (is.null(opt$output)){ print_usage(spec) }
if (is.null(opt$hb)){ print_usage(spec) }
if (is.null(opt$lb)){ print_usage(spec) }
if (is.null(opt$bs)){ opt$bs=30 }
print(opt$bs)
if (is.null(opt$ws)){ opt$ws=1e6 }

print(opt$ws)
if (is.null(opt$pta)){ opt$pta=95 }
print(opt$pta)
if (is.null(opt$ptb)){ opt$ptb=99 }
print(opt$ptb)
if (is.null(opt$pop)){ opt$pop="F2" }
print(opt$pop)

times<-Sys.time()
library(QTLseqr)

##ONE WAY
#library(magrittr)
#pwd<-setwd(opt$output)
#source(sprintf("%s%s",pwd,"/importFromGATK.R"))
#source(sprintf("%s%s",pwd,"/importfilter.R")
#source("/mnt/ilustre/users/meng.luo/QTLseqR/BSA_G_ED/bin/importFromGATK.R")
#source("/mnt/ilustre/users/meng.luo/QTLseqR/BSA_G_ED/bin/importfilter.R")
##TWO WAY

setwd(opt$output)

#chr<-opt$chr
#Chroms <- paste0(rep("Chr",chr), 1:chr)

#windowSize<- opt$ws

##Import SNP data from GATKTable file
HighBulk <- opt$hb
LowBulk <- opt$lb
data <-importFromGATK(file =opt$input,highBulk = opt$hb,lowBulk = opt$lb)

## if else for another way to filter data
#data<-filterSNPs(SNPset = data,refAlleleFreq = 0.20,minTotalDepth = 100,maxTotalDepth = 400,minSampleDepth = 40,minGQ = 99)

#write.table(data,file = "filter.index.result",sep = "\t",quote=F,row.names = FALSE)

seq_G <- runGprimeAnalysis(SNPset = data,windowSize =opt$ws,outlierFilter = "deltaSNP")

write.table(seq_G,file = "G.analysis.result",sep = "\t",quote=F,row.names = FALSE)

seq <- runQTLseqAnalysis(SNPset =data,windowSize = opt$ws,popStruc =paste(opt$pop),bulkSize =opt$bs,replications = 10000,intervals = c(opt$pta, opt$ptb))

write.table(seq,file = "qtlseq.analysis.result",sep = "\t",quote=F,row.names = FALSE)

library(parallel)
loessopt <- function (s, e, p, m) {
  x <- try(loess(e ~ p, span=s, degree=1, family="symmetric", surface='direct'), silent=T)
  if(class(x)=="try-error"){return(NA)}
  span <- x$pars$span
  n <- x$n
  traceL <- x$trace.hat
  sigma2 <- sum( x$residuals^2 ) / (n-1)
  delta1 <- x$one.delta
  delta2 <- x$two.delta
  enp <- x$enp
  if(m=="aicc"){return(log(sigma2) + 1 + 2* (2*(traceL+1)) / (n-traceL-2))}
  if(m=="gcv"){return(n*sigma2 / (n-traceL)^2)}
}
## three node
c1 <- makeCluster(as.numeric(3))
A<-seq_G$AD_REF.LOW/(seq_G$AD_REF.LOW+seq_G$AD_REF.HIGH)
B<-seq_G$AD_REF.HIGH/(seq_G$AD_REF.LOW+seq_G$AD_REF.HIGH)
C<-seq_G$AD_ALT.LOW/(seq_G$AD_ALT.HIGH+seq_G$AD_ALT.LOW)
D<-seq_G$AD_ALT.HIGH/(seq_G$AD_ALT.HIGH+seq_G$AD_ALT.LOW)
ED_AD<-sqrt((B-A)^2+(D-C)^2)
ED_DP<-sqrt((seq_G$DP.HIGH -seq_G$DP.LOW)^2)
mat<-as.matrix(cbind(seq_G$AD_REF.LOW,seq_G$AD_REF.HIGH,
                     seq_G$AD_ALT.LOW,seq_G$AD_ALT.HIGH))
n<-dim(mat)[1];m<-dim(mat)[2];ED_ADG<-NULL
for (i in 1:n) {
  ED_ADG[i]<-sqrt((mat[i,2]-mat[i,1]) ^ 2+(mat[i,4]-mat[i,3])^2)
}
ED_INDEX<-sqrt((seq_G$SNPindex.LOW-seq_G$SNPindex.HIGH)^2)
chrs <- unique(seq_G$CHROM)
#chrs <- chrs[grep(chrname, chrs)]
power<-4
m<-"aicc"
#ED<-ED_ADG
ED<-ED_INDEX
#ED<-ED_AD
#ED<-ED_DP
times<-Sys.time()
for(chr in chrs) {
  e<-ED[seq_G$CHROM==chr]^power
  p<-seq_G$POS[seq_G$CHROM==chr]
  spanresults.span <- NULL
  spanresults.aicc <- NULL
  if(as.integer(length(p)) < 100){
    cat(chr, " has less than 100 snps (n=", length(p), ")\n", sep="")
    next
  }
  ###change AIC for more running times
  spanresults.span <- seq(round(50/length(p),digits=3),2*round(50/length(p),digits=3), .001)
  spanresults.aicc <- parSapply(c1, spanresults.span, loessopt, e=e, p=p, m=m)
  usespan <- spanresults.span[spanresults.aicc==min(spanresults.aicc, na.rm=T)]
  lo <- loess(e ~ p, span = usespan[1], degree=1, family='symmetric', surface='direct') #designate usespan index incase of tie.
  seq_G$fitted[seq_G$CHROM==chr] <- lo$fitted
  #seq_G$unfitted[seq_G$CHROM==chr] <- lo$y
  #cat(chr, length(p), round(usespan[1], digits=3), "\n", sep="\t")
}
losefittime<-Sys.time()-times;print(losefittime);
seq_G$fitted[seq_G$fitted<0]<-0
cutoff <- 3*(sd(seq_G$fitted)+median(seq_G$fitted))
cutoff<-c(matrix(cutoff,length(seq_G$fitted)))
#seq_G$cutoff<-cutoff
ED4<-ED_AD^4
result<-cbind(seq_G[,c(1:4)],seq_G[,c(18,24,26,27,28)],seq[,c(31)],seq_G[,c(30,31,33)],ED4,seq_G[,c(34)],cutoff)
colnames(result)<-c("Chr","Pos","Ref","Alt","index.low",
	"index.high","delta","nSNP","WindowDeltaSNP",paste("I",95,sep=""),"Gprime","Gpvalue",
	"Gqvalue","ED4","ED4_fit","ED4_fit_cutoff")

write.table(result,file="three.methods.result",sep="\t",quote=F,row.names=F,col.names=T)

##plot






SNPset<-result

##ED

dat1<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$ED4_fit,SNPset$ED4_fit_cutoff)
colnames(dat1)<-c("Chr","Pos","ED4_fit","cutoff")
write.table(dat1,file="ED.result",sep="\t",quote=F,row.names=F,col.names=T)

##Gprime

dat2<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$Gprime,SNPset$Gqvalue)
colnames(dat2)<-c("Chr","Pos","Gprime","Gqvalue")
write.table(dat2,file="Gprime.result",sep="\t",quote=F,row.names=F,col.names=T)

##DeltaSNP
SNPset[,9]<-abs(SNPset[,9])
SNPset[,10]<-abs(SNPset[,10])

dat3<-data.frame(SNPset$Chr,SNPset$Pos,SNPset$WindowDeltaSNP,SNPset[,10])
colnames(dat3)<-c("Chr","Pos","WindowDeltaSNP",paste("I",95,sep=""))
write.table(dat3,file="Index.result",sep="\t",quote=F,row.names=F,col.names=T)



escaptime=Sys.time()-times;
print("Done!");
print(escaptime)

