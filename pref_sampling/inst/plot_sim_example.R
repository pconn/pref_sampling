### plot preferential sampling plots using data from "run_sim_example.R"
library(sp)
library(maptools)
library(rgeos)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)

setwd('c:/users/paul.conn/git/pref_sampling')
source('../SpatPred/SpatPred/R/util_funcs.R')
load('./Sim_data/sim_plot_data.Rdata')

N.true=matrix(Sim$N.true,ncol=1)
N.1=matrix(Sim$N.1,ncol=1)
N.2=matrix(Sim$N.2,ncol=1)
N.3=matrix(Sim$N.3,ncol=1)
C.1=matrix(Sim$C.1,ncol=1)
C.2=matrix(Sim$C.2,ncol=1)
C.3=matrix(Sim$C.3,ncol=1)
Delta=matrix(Sim$Delta,ncol=1)
Cov=matrix(Sim$Cov,ncol=1)

greenPalette <- colorRampPalette(brewer.pal(9, "Greens"))
bluePalette <- colorRampPalette(brewer.pal(9, "Blues"))
YlOrBr<- colorRampPalette(brewer.pal(9, "YlOrBr"))
txt.size=8


#turn sim data into SpatialPolygonsDataFrame for use with plot_N_map
Grid.topo=GridTopology(c(0,0),c(1,1),c(25,25))
Grid.SpG=SpatialGrid(Grid.topo)
Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(Grid.SpP)=CRS(laea_180_proj)   
Grid=SpatialPolygonsDataFrame(Grid.SpP,data=data.frame(Cov=Cov),match.ID=FALSE)



Grid.list=vector("list",1)
Grid.list[[1]]=Grid
max.N=100
Ntrue.plot=plot_N_map(1,N.true,Grid=Grid.list)+ggtitle(expression(paste("C. True abundance (",N[i],")")))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.1=plot_N_map(1,N.1,Grid=Grid.list)+ggtitle(expression(paste("G. ",hat(N[i]),', b=0')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.2=plot_N_map(1,N.2,Grid=Grid.list)+ggtitle(expression(paste("H. ",hat(N[i]),', b=1')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.3=plot_N_map(1,N.3,Grid=Grid.list)+ggtitle(expression(paste("I. ",hat(N[i]),', b=5')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
max.N=50
C.plot.1=plot_N_map(1,C.1,Grid=Grid.list)+ggtitle(expression(paste('D. Counts (',Y[i],'), b=0')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
C.plot.2=plot_N_map(1,C.2,Grid=Grid.list)+ggtitle(expression(paste('E. Counts (',Y[i],'), b=1')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
C.plot.3=plot_N_map(1,C.3,Grid=Grid.list)+ggtitle(expression(paste('F. Counts (',Y[i],'), b=5')))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
min.N=-2.5
max.N=2.5
Cov.plot=plot_N_map(1,Cov,Grid=Grid.list)+ggtitle("A. Covariate")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=greenPalette(100),limits=c(min.N,max.N))   
Delta.plot=plot_N_map(1,Delta,Grid=Grid.list)+ggtitle(expression(paste("B. Spatial random effects (",delta[i],")")))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=bluePalette(100),limits=c(min.N,max.N))   



pdf(file="Pref_samp_sim_maps.pdf")
grid.arrange(arrangeGrob(Cov.plot,Delta.plot,Ntrue.plot,C.plot.1,C.plot.2,C.plot.3,N.plot.1,N.plot.2,N.plot.3)) #,widths=unit(0.33,"npc"),heights=unit(0.33,"npc"),nrow=3))
dev.off()




    