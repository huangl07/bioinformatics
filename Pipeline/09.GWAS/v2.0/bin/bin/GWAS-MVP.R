library('getopt');
options(bitmapType='cairo')

spec = matrix(c(
	'hapmap','h',0,'character',
	'trait','t',0,'character',
	'output','o',0,'character',
	'help','m',0,'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("	
Usage:
	--hapmap	the input hapmap file
	--trait	the trait file 
	--output	the output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$trait)) { print_usage(spec)}
if ( is.null(opt$hapmap)){ print_usage(spec) }
if ( is.null(opt$output)){ print_usage(spec) }
times<-Sys.time()

library(MVP)
if(!dir.exists(opt$output)){dir.create(opt$output)}
setwd(opt$output)
MVP.Data(
	fileHMP="opt$hapmap",
	filePhe="opt$trait",
	sep.hmp="\t",
	sep.phe="\t",
	SNP.effect="Add",
	fileKin=FALSE, 
	filePC=FALSE,
	out="pop",
	priority="memory"
)
genotype <- attach.big.matrix("pop.geno.desc")
phenotype <- read.table("pop.phe",head=TRUE)
map <- read.table("pop.map" , head = TRUE)
for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    nPC.GLM=5,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    perc=1,
    priority="memory",
    ncpus=16,
    vc.method="EMMA",
    maxLoop=10,
    method.bin="FaST-LMM",
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU"),
	file="pdf",
  )
}

escaptime=Sys.time()-times;
print("Done!")
print(escaptime)
