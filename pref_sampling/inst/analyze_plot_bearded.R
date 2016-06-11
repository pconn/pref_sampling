### analyze and visualize bearded seal model output
library(ggplot2)
library(grid)
library(gridExtra)
setwd('c:/users/paul.conn/git/pref_sampling')
load("AlaskaBeringData2012_2013_14Dec2015.Rdat")  #from bass package on github/nmml
source('./pref_sampling/R/util_funcs.R')
setwd('c:/users/paul.conn/git/pref_sampling/Output')
n.cells=nrow(Data$Adj)
  
N.table = matrix(0,6,4)
colnames(N.table)=c("cov.mod","b.mod","N.est","SE")
AIC.table = matrix(0,6,5)
colnames(AIC.table)=c("cov.mod","b.mod","Likelihood","No.params","AICc")

N.b.par=c(0,1,3)
N.cov.par=c(0,3)

N.s=matrix(0,n.cells,6)

n.counter=1
Plot=vector("list",1)
for(icov in 1:2){
  for(ib in 1:3){
    fname=paste0("OutBearded_cov",icov,"_b",ib,".RData")
    load(file=fname)
    N.table[n.counter,1:2]=AIC.table[n.counter,1:2]=c(icov,ib)
    N.table[n.counter,3]=Out$SD$unbiased$value["total_abundance"]
    cur.pl=which(names(Out$SD$unbiased$value)=="total_abundance")
    N.table[n.counter,4]=Out$SD$sd[cur.pl]
    AIC.table[n.counter,3]=Out$log_lik
    AIC.table[n.counter,4]=N.b.par[ib]+N.cov.par[icov]*2+3
    AIC.table[n.counter,5]=-2*AIC.table[n.counter,3]+2*AIC.table[n.counter,4]
    
    N.s[,n.counter]=Out$SD$unbiased$value[which(names(Out$SD$unbiased$value)=="Z_s")]
    Plot[[n.counter]]=plot_N_map(cur.t=1,N=matrix(N.s[,n.counter],ncol=1),Grid=Data$Grid$y2012[[1]])

    n.counter=n.counter+1
  }
}

#Make some plots to examine response surface for preferential sampling
Grid.new=Data$Grid$y2012[[1]]
n_cells=length(Grid.new)

fname=paste0("OutBearded_cov",2,"_b",2,".RData") #constant b(s)
load(file=fname)
Delta=Out$SD$par.random[1:1331]
Eta=Out$Report$eta_s[1:1331]
b=Out$Report$beta_b
Eta22_plot = plot_N_map(cur.t=1,N=matrix(Eta,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="eta(s)")
Delta22_plot = plot_N_map(cur.t=1,N=matrix(Delta,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="delta(s)")
bDelta22_plot = plot_N_map(cur.t=1,N=matrix(Delta*b,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="b*delta(s)")

fname=paste0("OutBearded_cov",2,"_b",3,".RData") #constant b(s)
load(file=fname)
x_s_b = matrix(1,nrow=n_cells,ncol=3)  #design matrix for preferential sampling effect
Centroids=gCentroid(Grid.new,byid=TRUE)
loc_s=cbind(Centroids$x,Centroids$y)/25067.53 #standardize by making distance between neighboring cells = 1.0
colnames(loc_s)=c("x","y")
rownames(loc_s)=c(1:n_cells)
x_s_b[,2:3] = loc_s
Pref_effect=x_s_b%*%matrix(Out$SD$par.fixed[8:10],ncol=1)
Delta=Out$SD$par.random[1:1331]
Eta=Out$Report$eta_s[1:1331]
Eta23_plot = plot_N_map(cur.t=1,N=matrix(Eta,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="eta(s)")
Delta23_plot = plot_N_map(cur.t=1,N=matrix(Delta,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="delta(s)")
bDelta23_plot=plot_N_map(cur.t=1,N=matrix(Delta*Pref_effect,ncol=1),Grid=Data$Grid$y2012[[1]],leg.title="b(s)*delta(s)")
Eta22_plot=Eta22_plot+ggtitle("Cov=1, b=1")
Eta23_plot=Eta23_plot+ggtitle("Cov=1, b=2")
pdf("bearded_RE_maps.pdf")
grid.arrange(arrangeGrob(Eta22_plot,Eta23_plot,Delta22_plot,Delta23_plot,bDelta22_plot,bDelta23_plot,widths=rep(unit(0.5,"npc"),2),heights=rep(unit(0.33,"npc"),3),nrow=3))
dev.off()


max.N=max(N.s)
library(RColorBrewer)
myPalette=colorRampPalette(brewer.pal(9, "Purples")) 
for(i in 1:6)Plot[[i]]=Plot[[i]]+scale_fill_gradientn(colours=myPalette(100),limits=c(0,max.N))

Plot[[1]]=Plot[[1]]+ggtitle("A. Cov = 0; b = 0")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+ylab('Northing')+xlab('')
Plot[[4]]=Plot[[4]]+ggtitle("B. Cov = 1; b = 0")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')
Plot[[2]]=Plot[[2]]+ggtitle("C. Cov = 0; b = 1")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+ylab('Northing')+xlab('')
Plot[[5]]=Plot[[5]]+ggtitle("D. Cov = 1; b = 1")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')
Plot[[3]]=Plot[[3]]+ggtitle("E. Cov = 0; b = 2")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+ylab('Northing')+xlab('Easting')
Plot[[6]]=Plot[[6]]+ggtitle("F. Cov = 1; b = 2")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('Easting')+ylab('')

pdf("bearded_maps.pdf")
grid.arrange(arrangeGrob(Plot[[1]],Plot[[4]],Plot[[2]],Plot[[5]],Plot[[3]],Plot[[6]],widths=rep(unit(0.5,"npc"),2),heights=rep(unit(0.33,"npc"),3),nrow=3))
dev.off()

