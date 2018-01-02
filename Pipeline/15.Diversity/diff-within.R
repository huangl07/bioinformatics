library('getopt');
options(bitmapType='cairo')
spec = matrix(c(
	'input','v',0,'character',
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
	--input	the genalex file
	--out	the output file
	--help		usage
\n")
	q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$input)){ print_usage(spec) }
if ( is.null(opt$out)){ print_usage(spec) }

times<-Sys.time()
library(poppr)
library(ggplot2)
library(reshape2)

monpop <- read.genalex(opt$input)
setPop(monpop)<- ~Pop
monpop<-missingno(monpop, type = "mean", cutoff = 0.05, quiet = FALSE, freq = FALSE)
popdiversity<-poppr(monpop)
write.table(popdiversity,file=opt$out,sep="\t",row.names=FALSE)

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
