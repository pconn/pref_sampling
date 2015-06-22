#' function to simulate abundance over a simulated landscape
#' @param S Number of cells (must be a square number)
#' @param n.covs Number of habitat covariates (default 4)
#' @param tau.epsilon Precision of exchangeable errors in log space
#' @return A SpatialPolygonsDataFrame holding covariates, abundance, and lambda (expected abundance)
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic<-function(S,n.covs=4,tau.epsilon=100){
  DEBUG=FALSE   
  PLOT=FALSE
  if(PLOT)set.seed(11111)
  tau.epsilon=100
  #if(DEBUG)tau.epsilon=100
  
  #source('./SpatPred/R/util_funcs.R')
   
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
  
  x.len=sqrt(S)
 
  n.basis=10
  n.zero=6
  Omega=matrix(0,n.basis,S)
  Cov=matrix(0,n.covs,S)
  flag=1
  while(flag==1){  #reject if any |cor|>0.75
    for(i in 1:n.basis){
      my.mod=RMexp(var=1,scale=runif(1,5,100))
      Omega[i,]=RFsimulate(model=my.mod,x=c(1:x.len),y=c(1:x.len))$variable1
    }
    for(i in 1:n.covs){
      Wts=matrix(rnorm(n.basis),1,n.basis)
      Wts[sample(c(1:n.basis),n.zero)]=0
      #cat(Wts)
      #cat('\n')
      Cov[i,]=Wts%*%Omega
    }
    Covs=t((Cov-rowMeans(Cov))/sqrt(apply(Cov,1,'var')))
    cur.cor=max(abs(c(cor(Covs[,1],Covs[,2]),cor(Covs[,1],Covs[,3]),cor(Covs[,2],Covs[,3]))))
    if(cur.cor<0.75)flag=0
  }
    
  
  if(PLOT==TRUE){
    #XY=data.frame(x=Dat.matern$x,y=Dat.matern$y)
    Grid1=data.frame(y=rep(1:30,each=30),x=rep(c(1:30),30),cov=as.vector(Covs[,1]))
    Grid2=data.frame(y=rep(1:30,each=30),x=rep(c(1:30),30),cov=as.vector(Covs[,2]))
    Grid3=data.frame(y=rep(1:30,each=30),x=rep(c(1:30),30),cov=as.vector(Covs[,3]))
    Grid=rbind(Grid1,Grid2,Grid3)
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
     Grid.plot1=ggplot()+geom_raster(data=Grid1,aes(x=x,y=y,fill=cov))+theme(title=element_text(size=rel(1.5)),axis.text=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0))+ggtitle("A.")+scale_fill_gradientn(colours=myPalette(100))
    Grid.plot2=ggplot()+geom_raster(data=Grid2,aes(x=x,y=y,fill=cov))+theme(title=element_text(size=rel(1.5)),axis.text=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0))+ggtitle("B.")+scale_fill_gradientn(colours=myPalette(100))
    Grid.plot3=ggplot()+geom_raster(data=Grid3,aes(x=x,y=y,fill=cov))+theme(title=element_text(size=rel(1.5)),axis.text=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0))+ggtitle("C.")+scale_fill_gradientn(colours=myPalette(100))    
    #now, iterate to produce a histogram of implied correlation coefficients
    Cor=rep(0,1000)
    for(irep in 1:1000){
      my.mod=RMexp(var=1,scale=100)
      for(i in 1:n.basis){
        my.mod=RMexp(var=1,scale=runif(1,5,100))
        Omega[i,]=RFsimulate(model=my.mod,x=c(1:x.len),y=c(1:x.len))$variable1
      }
      for(i in 1:2){
        Wts=matrix(rnorm(n.basis),1,n.basis)
        Wts[sample(c(1:n.basis),n.zero)]=0
        Cov[i,]=Wts%*%Omega
      }
      Covs=t(Cov/apply(Cov,1,'max'))
      Cor[irep]=cor(Covs[,1],Covs[,2])
    }
    Cor.df=data.frame(Cor=Cor)
    
    Cov.plot=ggplot()+geom_density(data=Cor.df,aes(x=Cor))+xlab("Correlation coefficient")+ylab("Density")+
      theme(plot.title = element_text(hjust = 0),axis.text=element_text(size=rel(1.5)),title=element_text(size=rel(1.5)))+ggtitle("D.")
    pdf(file="Covs_coregion.pdf")
    library(gridExtra)
    grid.arrange(arrangeGrob(Grid.plot1,Grid.plot2,Grid.plot3,Cov.plot,nrow=2))
    dev.off()
  }

  col.names="cov1"
  if(n.covs>1){
    for(i in 2:n.covs)col.names=c(col.names,paste("cov",i,sep=''))
  }
  colnames(Covs)=col.names
  Covs=as.data.frame(Covs)
  
  
  Grid.topo=GridTopology(c(0,0),c(1,1),c(x.len,x.len))
  Grid.SpG=SpatialGrid(Grid.topo)
  Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
  Grid.SpPDF=SpatialPolygonsDataFrame(Grid.SpP,data=Covs,match.ID=FALSE)
  laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                         "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(Grid.SpPDF)=CRS(laea_180_proj)   
  
  N.sum=100000
  while(N.sum>99999){
    Beta=c(rnorm(1,2.5,0.5),rnorm(n.covs,0,0.4),rnorm(n.covs,0,0.2),rnorm(choose(n.covs,2),0,0.2))
    X=model.matrix(~poly(as.matrix(Covs),degree=2,raw=TRUE))
    Lambda=exp(X%*%Beta+rnorm(S,0,sqrt(1/tau.epsilon)))
    lam.20=sum(Lambda[which(rank(Lambda)>880)])
    if(is.na(sum(Lambda))==FALSE & max(Lambda)<100000 & lam.20/sum(Lambda)<0.9){
      N=rpois(S,Lambda)
      N.sum=sum(N)
    }
  }
  
  Data=Grid.SpPDF
  Data@data=data.frame(Easting=rep(1:sqrt(S),sqrt(S)),Northing=rep(1:sqrt(S),each=sqrt(S)),Lambda=Lambda,N=N)
  Data@data=cbind(Data@data,as.data.frame(Covs))
  
  
  #p1=ggplot(Data)+aes(Easting,Northing,fill=N)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov1)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov2)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov3)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov4)+geom_raster()
  #p1   
  return(Data)
}


#' function to simulate spatial count data over a simulated grid, conditional on expected abundance
#' @param Data A SpatialPolygonsDataFrame for a square survey grid including a data column "N" that gives total abundance in each survey unit
#' @param S Number of cells (must be a square number)
#' @param type Type of spatial survey.  Choices are "random" (default), "clustered" (at center of grid), "balanced" (spatially balanced using random tesselation design), or "preferential" (concentrating on areas of high abundance)
#' @param prop.sampled A numeric value giving the proportion of each survey unit that is actually sampled (default 0.1)
#' @param n.points Number of survey units that are sampled.  The default is to round (0.1*S)
#' @param pref.mult For type="preferential", this determines how much more likely it is to sample the highest abundance than the lowest abundance
#' @return A data.frame with two columns, "Mapping" - indicating which survey units are sampled, and "Counts" - giving no. of animals counted
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_effort<-function(Data,S,type="random",prop.sampled=0.1,n.points=round(0.1*S),pref.mult=2){
  x.len=sqrt(S)
  if(type=="preferential"){
    Rank=rank(Data[["N"]])
    Inc.prob=1+c(1:S)*(pref.mult-1)/(S-pref.mult)
    Inc.prob=Inc.prob[Rank]
    Mapping=sample(c(1:S),n.points,prob=Inc.prob)
  }
    
  if(type=="clustered"){
    center=sqrt(S)/2
    XY=expand.grid(x=c(1:sqrt(S)),y=c(1:sqrt(S)))
    Distances=sqrt((XY[,1]-center)^2+(XY[,2]-center)^2)
    Inc.prob = exp(-0.2*Distances)
    Inc.prob = Inc.prob/sum(Inc.prob)
    Mapping=sample(c(1:S),n.points,prob=Inc.prob)

    #ids=factor(c(1:n.points))
    #df<-data.frame(id=ids,x=XY[Mapping,"x"],y=XY[Mapping,"y"])
    #Abund=runif(S,0,100) #sum=11971; expected total abundance =11971*4=47885
    #Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
    #colnames(Abund.df)=c("y","x","Abundance")
    #require(ggplot2)
    #plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
    #plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
    #plot1    
  }
  if(type=="balanced"){
    design<-list("Stratum1"=list(panel=c(Panel=n.points),seltype="Equal",over=0))
    att.frame=data.frame(x=rep(c(1:sqrt(S)),each=sqrt(S)),y=rep(c(1:(sqrt(S))),sqrt(S)))
    spsurv<-grts(design=design,src.frame="att.frame",att.frame=att.frame,shapefile=FALSE)
    
    x.base=spsurv$x
    y.base=spsurv$y #c(spsurv$y*2-1,spsurv$y*2)
    ids=factor(c(1:n.points))
    df<-data.frame(id=ids,x=x.base,y=y.base)
    
    #plot expected abundance
    #Abund=4*Exp.grp.abund #sum=11971; expected total abundance =11971*4=47885
    #Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
    #colnames(Abund.df)=c("y","x","Abundance")
    #require(ggplot2)
    #plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
    #plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
    #plot1
    
    #Determine mapping, fraction of cell occupied by each transect
    Mapping=rep(0,n.points)
    for(itrans in 1:n.points){
      Mapping[itrans]=(x.base[itrans]-1)*sqrt(S)+sqrt(S)-y.base[itrans]+1
    }    
  }  
  if(type=="random"){
    Mapping=sort(sample(c(1:S),n.points,replace=FALSE))
  }
  Counts=rbinom(n.points,Data@data[,"N"][Mapping],prop.sampled)
  Effort=data.frame(Mapping=Mapping,Counts=Counts)  
}  
