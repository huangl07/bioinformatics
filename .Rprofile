local({r <- getOption("repos")
	r["CRAN"] <- "http://mirrors.ustc.edu.cn/CRAN/"
	r["CRANextra"] <- "http://mirrors.aliyun.com/CRAN/"
	options(repos=r)})
bioclite<-function(package)
{
	source("https://bioconductor.org/biocLite.R")
	options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
	biocLite(package)
}
