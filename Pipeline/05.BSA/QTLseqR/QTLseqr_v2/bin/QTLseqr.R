#!/usr/bin/env Rscript
times<-Sys.time()
library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','i',1,'character',
	'output','o',0,'character',
	'HB','p',1,'character',
	'LB','l',1,'character',
	'help','h',0,'logical'
	), byrow=TRUE, ncol=4)
opt = getopt(spec)
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage example: 
	Rscript emap.R --input --output --HB --LB
Usage:
	--input *vcf.table
	--output	output dir
	--HB	HighBulk sampleID
	--LB	LowBulk sampleID
	--help		usage
\n")
	q(status=1);
}

times<-Sys.time()


if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)) { print_usage(spec)}
if ( is.null(opt$output)){ print_usage(spec) }
if ( is.null(opt$HB)){ print_usage(spec) }
if ( is.null(opt$LB)){ print_usage(spec) }

library("QTLseqr")
setwd(opt$output)
HighBulk = opt$HB
LowBulk = opt$LB

df <-importFromGATK(file = opt$input,highBulk = HighBulk,lowBulk = LowBulk)

df_filt <-filterSNPs(SNPset = df,refAlleleFreq = 0.20,minTotalDepth = 10,maxTotalDepth = 100,minSampleDepth = 10,minGQ = 99)
write.table(df_filt,file = "filter.index.result",sep = "\t",row.names = FALSE)
	
df_filt <- runGprimeAnalysis(SNPset = df_filt,windowSize = 1e6,outlierFilter = "deltaSNP")
write.table(df_filt,file = "G.analysis.result",sep = "\t",row.names = FALSE)

df_filt <- runQTLseqAnalysis(SNPset = df_filt,windowSize = 1e6,popStruc = "F2",bulkSize = 25,replications = 10000,intervals = c(95, 99))
write.table(df_filt,file = "qtlseq.analysis.result",sep = "\t",row.names = FALSE)

png(file = "bsa.Gprime.png")
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.001)
dev.off()
pdf(file = "bsa.Gprime.pdf")
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.001)
dev.off()
png(file = "bsa.deltaSNP.png")
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()
pdf(file = "bsa.deltaSNP.pdf")
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()
png(file = "bsa.negLog10Pval.png")
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotIntervals = TRUE, q = 0.001)
dev.off()
pdf(file = "bsa.negLog10Pval.pdf")
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotIntervals = TRUE, q = 0.001)
dev.off()
png(file = "bsa.nSNPs.png")
plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)
dev.off()
pdf(file = "bsa.nSNPs.pdf")
plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)
dev.off()
getQTLTable(SNPset = df_filt, alpha = 0.001, export = TRUE, fileName = "qtl.csv")