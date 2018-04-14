######################################## creat Manhattan plot

rm(list=ls())
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")
library("EMMREML")
library(ape)
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
    ha2=c(0.8,0.8)
    rg=diag(1,2)
    rg[1,2]=rg[2,1]=0.9
    re=diag(1,2)
    re[1,2]=re[2,1]=.1
    NQTN=20
    NE=2
rep=50
n_e=NE
set.seed(99163)

system(paste("./GEMT --gwas --file ","mygemtData --numeric ","--interaction ",n_e,sep=""))
result=read.table(paste("GEMT","_GWAS_result.txt",sep=""),head=T)
result=result[,c(1,2,3,5,4)]
chrom=result[1:(nrow(result)/n_e),2]
for(i in 1:n_e)
{
	result[(1+(i-1)*(nrow(result)/n_e)):((nrow(result)/n_e+(i-1)*(nrow(result)/n_e))),2]=chrom
}


source("/Users/Jiabo/Documents/R/code/Manhattan_g&e.R")
#result=cbind(newmap,result[,4])
GAPIT.Manhattan(GI.MP=result[,-1],name.of.trait="GW",Name_environ=c("13_FC_Irr","13_FC_Dry","12_Gr_Irr","12_Gr_Dry"),Nenviron=n_e)



######################################## complare power with blink_mean MLM_mean

rm(list=ls())
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")
library("EMMREML")
library(ape)
source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
    ha2=c(0.8,0.8)
    rg=diag(1,2)
    rg[1,2]=rg[2,1]=0.9
    re=diag(1,2)
    re[1,2]=re[2,1]=.1
    NQTN=20
    NE=2
rep=50
n_e=NE
set.seed(99163)
###########################################
#setwd("/Users/Jiabo/Documents/GEMT/281demo")

source("/Users/Jiabo/Documents/R/code/G_E.Duplicated.R")
source("/Users/Jiabo/Documents/R/code/G&E_Simulation.R")
source("/Users/Jiabo/Dropbox/GAPIT/Functions/GAPIT.FDR.TypeI.R")
source("/Users/Jiabo/Dropbox/GAPIT/Functions/GAPIT.AUC.R")

GD=read.table("gemtData.dat",head=F)
GM=read.table("gemtData.map",head=T)

taxa=NULL
E_name=NULL
store_map=GM
for(i in 1:ncol(GD))
{
taxa[i]=paste("taxa_",i,sep="")

}

mean_fdr=NULL
dupl_fdr=NULL
gemt_fdr=NULL
add_fdr=NULL
gemt_fdr2=NULL

for(k in 1:rep)
{

mysimu=G_E.Simulation(GD=cbind(as.data.frame(taxa),t(GD)),ha2=ha2,NQTN=NQTN,NE=NE,rg=rg,re=re)
myY=mysimu$Y
cor(myY[,2],myY[,3])


#addposition=mysimu$QTN.position[c(1:NQTN)]
#interposition=mysimu$QTN.position[-c(1:NQTN)]
allposition=mysimu$QTN.position
Y=apply(as.matrix(myY[,-1]),1,mean)
Y=cbind(as.data.frame(taxa),Y)
colnames(Y)=c("taxa","blink")
dupli_taxa=NULL
dupli_Y=NULL
GDGD=GD
mapmap=NULL
#dpuliGD=GD
#dpuliGD[]=0
store_map=GM

diag_matrix=matrix(c(1,1,1,0),2,2)
#diag_matrix[1,]=1


X=as.matrix(GD)
dupliGD=kronecker(diag_matrix,X)

for(i in 1:NE)
{
	dupli_taxa=append(dupli_taxa,taxa)
    dupli_Y=append(dupli_Y,myY[,i+1])
    if(i!=NE)GDGD=cbind(GDGD,GD)

    store_map[,2]=store_map[,2]+max(store_map[,2])*(i-1)
    mapmap=rbind(mapmap,store_map)
}
YY=cbind(as.data.frame(dupli_taxa),dupli_Y)
colnames(YY)=c("taxa","dupli_trait")
#for GEMT
write.table(myY,"mygemtData.txt",,row.names=F,quote=F)
write.table(GD,"mygemtData.dat",,row.names=F,quote=F,col.names=F)
write.table(GM,"mygemtData.map",,row.names=F,quote=F,col.names=T)

#for blink duplicated
write.table(YY,"myblink_dupData.txt",,row.names=F,quote=F)
#write.table(GD,"wheat_10K.dat",,row.names=F,quote=F,col.names=F)
write.table(GDGD,"myblink_dupData.dat",,row.names=F,quote=F,col.names=F)
write.table(GM,"myblink_dupData.map",,row.names=F,quote=F,col.names=T)
#for blink mean
write.table(Y,"blinkData.txt",,row.names=F,quote=F)

system("cp mygemtData.dat blinkData.dat")
system("cp mygemtData.map blinkData.map")

system(paste("./GEMT --gwas --file ","mygemtData --numeric ","--interaction ",n_e,sep=""))
system(paste("./biogpu0.9 --gwas --file ","myblink_dupData --numeric ",sep=""))
system(paste("./biogpu0.9 --gwas --file ","blinkData --numeric ",sep=""))
#system(paste("./biogpu0.9 --gwas --file ","gemtData --numeric ",sep=""))

N=nrow(GM)
mmGD=cbind(as.data.frame(Y[,1]),t(GD))
colnames(mmGD)=c("taxa",GM[,1])
mmGM=GM

#dpuliGD=GD
#for MLM
myadd <- GAPIT(
Y=Y,
GD=mmGD,
GM=mmGM,
#CV=myCV,
PCA.total=3,
group.from=100000,
group.to=100000,
file.output=F
)
aa=myadd$GWAS
bb=merge(mmGM,aa[,c(1,4,5)],x.by="SNP",y.by="SNP")
fdr_add=GAPIT.FDR.TypeI(GWAS=bb,GM=mmGM,WS=1,MaxBP=N,seqQTN=allposition)
add_fdr=cbind(add_fdr,fdr_add$FDR)

gemt_result=read.table("GEMT_GWAS_result.txt",head=T)
gemt_result=gemt_result[,c(1,2,3,5,4)]
fdr_gemt=GAPIT.FDR.TypeI(GWAS=gemt_result[1:nrow(GM),],GM=GM,WS=1,maxOut=N,seqQTN=allposition)
gemt_fdr=cbind(gemt_fdr,fdr_gemt$FDR)
fdr_gemt2=GAPIT.FDR.TypeI(GWAS=gemt_result[-c(1:nrow(GM)),],GM=GM,WS=1,maxOut=N,seqQTN=allposition)
gemt_fdr2=cbind(gemt_fdr2,fdr_gemt2$FDR)

blink_result=read.table("blink_GWAS_result.txt",head=T)
blink_result=blink_result[,c(1,2,3,5,4)]
fdr_mean=GAPIT.FDR.TypeI(GWAS=blink_result,GM=GM,WS=1,maxOut=N,seqQTN=allposition)
mean_fdr=cbind(mean_fdr,fdr_mean$FDR)



blink_result=read.table("dupli_trait_GWAS_result.txt",head=T)
blink_result=blink_result[,c(1,2,3,5,4)]
fdr_dupl=GAPIT.FDR.TypeI(GWAS=blink_result,GM=GM,WS=1,maxOut=N,seqQTN=allposition)
dupl_fdr=cbind(dupl_fdr,fdr_dupl$FDR)

}
#system("remove ")
gemt_plot=apply(gemt_fdr,1,mean)
mean_plot=apply(mean_fdr,1,mean)
dupl_plot=apply(dupl_fdr,1,mean)
add_plot=apply(add_fdr,1,mean)
gemt_plot2=apply(gemt_fdr2,1,mean)


add_power=seq(1/(NQTN),1,1/(NQTN))

tiff(file=paste("Mean_Duplicate_blink_GEMT_MLM_model_blink_",NE,"_Power against FDR.tif",sep=""),width=150, height=150, units='mm',res=200,compression='lzw')
par(mar=c(4,4,2.5,2.5))
plot(as.numeric(mean_plot),add_power,col="darkred",xlim=c(0,1),ylim=c(0,1),yaxp=c(0,1,10),xaxp=c(0,1,10),type="b",ylab="",xlab="",lwd=1.3)
par(new=T)
plot(as.numeric(dupl_plot),add_power,col="darkblue",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F,type="b",lwd=1.3)
par(new=T)
plot(as.numeric(gemt_plot),add_power,col="black",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F,type="b",lwd=1.3)
par(new=T)

plot(as.numeric(add_plot),add_power,col="yellow",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F,type="b",lwd=1.3)
par(new=T)
plot(as.numeric(gemt_plot2),add_power,col="darkgreen",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=F,type="b",lwd=1.3)


mtext("FDR",side=1,line=2.5,col="black")
mtext("POWER",side=2,line=2.5,col="black")
legend("topleft",legend=c("Mean_Blink","Duplicate_Blink","Mean_MLM","GEMT_add","GEMT_inter"),ncol=1,
col=c("darkred","darkblue","yellow","black","darkgreen"),pch=c(15,15,15),lty=0,lwd=1.3,cex=1,
bty = "n", bg = par("bg"))

dev.off()


