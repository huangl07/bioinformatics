library(qtl)
mapthis<-read.cross("csv",file="mapthis.csv",estimate.map=FALSE)
#miss
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < nind(mapthis)*0.7])
mapthis <- drop.markers(mapthis, todrop)
#seg
gt <- geno.table(mapthis)
todrop <- rownames(gt[gt$P.value < 0.05,])
mapthis <- drop.markers(mapthis, todrop)
#indi
cg <- comparegeno(mapthis)
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
mapthis <- subset(mapthis, ind=-wh[,2])
#group
mapthis <- est.rf(mapthis)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
#adjust 
toswitch<-markernames(mapthis,chr=c(1,9:11))
mapthis <- switchAlleles(mapthis, toswitch)
mapthis <- est.rf(mapthis)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
toswitch <- markernames(mapthis,chr=c(4,6))
mapthis <-switchAlleles(mapthis,toswitch)
mapthis <- est.rf(mapthis)
mapthis <- formLinkageGroups(mapthis, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
#mapping
mapthis <- orderMarkers(mapthis)
plotMap(mapthis)
plotRF(mapthis)
plotGeno(mapthis,chr=5)
