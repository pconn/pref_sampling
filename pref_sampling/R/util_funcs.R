 
#' Produce an RW1 adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj <- function(x,y,byrow=FALSE){
  Ind=matrix(c(1:(x*y)),y,x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,x*y,x*y)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i==1 & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i==1 & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i>1 & i<n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i>1 & i<n.row & j>1 & j<n.col){
        cur.nums=c(Ind[i,j]-n.row-1,Ind[i,j]-n.row,Ind[i,j]-n.row+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+n.row-1,Ind[i,j]+n.row,Ind[i,j]+n.row+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
        
      }
      if(i==n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
      }
      if(i==n.row & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
      if(i==n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
    }
  }
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

  
#' Function to plot abundance map for Bering Sea survey grid
#' @param cur.t  Time step to plot map for
#' @param N A matrix holding abundance point estimates; different rows correspond to different sampling units
#' @param Grid A spatial polygons data frame holding grid entries
#' @param highlight If provided, this vector specifies which cells to separately highlight during plotting
#' @param cell.width if highlight is provided, cell.width must provide the width of a composite grid cell
#' @param leg.title Title for legend of plot (if different from the default "Abundance")
#' @param hcolor If highight provided, gives color for highlighting (default = "yellow")
#' @param cur.max If provided, sets maximum value for plotted abundance
#' @param myPallete A colorRampPalette; default is a spectral palette from RColorBrewer
#' @return An abundance map
#' @export 
#' @keywords abundance map, plot
#' @author Paul Conn \email{paul.conn@@noaa.gov}
plot_N_map<-function(cur.t=1,N,Grid,highlight=NULL,cell.width=1,leg.title="Abundance",hcolor='yellow',cur.max=NULL,myPalette=NULL){
  #require(ggplot2)
  #require(rgeos)
  #require(RColorBrewer)
  #if(is.null(myPalette)==TRUE)myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  if(is.null(myPalette)==TRUE)myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  Tmp=Grid
  if(is.null(highlight)==FALSE){
    midpoints=data.frame(gCentroid(Tmp[highlight,],byid=TRUE))
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N
  Cur.df=cbind(data.frame(gCentroid(Tmp,byid=TRUE)),Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(is.null(cur.max))p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  else p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title,limits=c(0,cur.max))
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour=hcolor,aes(xmin=Easting-0.5*cell.width,xmax=Easting+0.5*cell.width,ymin=Northing-0.5*cell.width,ymax=Northing+0.5*cell.width))
  }
  p1
}

#' SIMULATE AN ICAR PROCESS 
#' @param Q Precision matrix for the ICAR process
#' @return Spatial random effects
#' @export 
#' @keywords ICAR, simulation
#' @author Devin Johnson
rrw <- function(Q){
  v <- eigen(Q, TRUE)
  val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
  P <- v$vectors
  sim <- P%*%diag(val.inv)%*%rnorm(dim(Q)[1], 0, 1)
  X <- rep(1,length(sim))
  if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))
  sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
  return(sim)
}

#' inverse logit transformation
#' @param x quantity to be transformed
#' @return real valued quantities
#' @export 
#' @keywords logit, expit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
expit<-function(x)1/(1+exp(-x))

#' logit transformation
#' @param x quantity to be transformed
#' @return logit value(s)
#' @export 
#' @keywords logit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
logit<-function(x)log(x/(1-x))


