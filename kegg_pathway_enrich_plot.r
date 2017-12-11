#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript

library(ggplot2)
# get args
args <-commandArgs(TRUE)



# reading data
ori_data <- read.delim(args[1], row.names = 1, header=TRUE)
head(ori_data)


# create data frame
# check pathway number
row_num <- min(dim(ori_data)[1], 20) 	# only drawing first 20 row
if(row_num<1){
	stop("the number of pathway < 1")
}
# create
df <- data.frame(rich_factor=ori_data[1:row_num, 2], fdr=-log10(ori_data[1:row_num, 3]))
KEGG_pathway <- factor(rownames(ori_data[1:row_num,]))


# plot
p <- qplot(df[,1], df[,2], colour=KEGG_pathway, shape=KEGG_pathway, size=I(4), xlab="enrichment factor", ylab="-log10(Q-value)") + scale_shape_manual(values=1:20)


# output plot
png(filename=args[2], height = 3000, width = 5000, res = 500, units = "px")
print(p)
dev.off()




