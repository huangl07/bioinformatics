times<-Sys.time() 

library('getopt')
options(bitmapType='cairo')
spec = matrix(c(
    'infile','f',0,'character',
    'outfile','o',0,'character',
    'order','l',0,'character'
     ), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
    cat(getopt(spec, usage=TRUE));
    cat("Usage example: \n")
    cat("
    Usage example: 

    Usage:
    --infile input file name key
    --order input order file
    --outfile output file name key 
    --help      usage
    \n")
    q(status=1);
}
if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$infile) ) { print_usage(spec) }
if ( is.null(opt$outfile) ) { print_usage(spec) }
source("/mnt/ilustre/users/nanshan.yang/newmdt/Pipline/03.history/treemix/bin/bin/plotting_funcs.R")

png(paste(opt$outfile,".png",sep=""),height=800*3,width=800*3,res=72*3)
plot_tree(opt$infile)
plot_resid(opt$infile,opt$order)
dev.off()
pdf(paste(opt$outfile,".pdf",sep=""),height=800,width=800)
plot_tree(opt$infile)
plot_resid(opt$infile,opt$order)
dev.off()
escaptime=Sys.time()-times;
print("Done!")
print(escaptime) 
