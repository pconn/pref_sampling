#### plot expectations of probability of sampling vs. abundance residuals


# Z -- true abundance
# Y -- sampled abundance
# R -- predicted probability of sampling
# delta -- affects density
# eta -- affects sampling probability

library( RandomFields )
library( INLA )

# change directory
TmbDir = "C:/Users/paul.conn/git/pref_sampling/TMB_version/inst/executables"
setwd( TmbDir )

# SOurce
source( paste0(TmbDir,"/../../examples/rect_adj.R") )
source( paste0(TmbDir,"/../../examples/rrw.R") )

# Settings
grid_dim = c("x"=25, "y"=25)
n_samp = 50
n_cells = grid_dim[1]*grid_dim[2]
SpatialScale = sqrt(prod(grid_dim))/5  # Range ~ 2*Scale
SD_eta = SD_x = SD_delta = 1
beta0 = 2
Alpha = 1  # Smoothness for GMRF, 1 or 2 (1 is faster)
n_sim = 500

# Configurations
b_set = c(0,1,5)            # Impact of delta on sampling intensity
Data=matrix(0,4,n_sim*n_cells)

# Loop through replicates
set.seed(12345)
irow=1
for(i in 1:n_sim){
  #simulate abundance
  # Settings

  #note: betax set to zero in last commit by Jim (v4)
  betax = runif(1,-.5,.5)       # Impact of x on density
  betax_prob = runif(1,-.5,.5)   # Impact of x on sampling intensity
  #betax=0
  #betax_prob=0
  
  # Spatial model
  loc_s = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
  model_delta <- RMgauss(var=SD_delta^2, scale=SpatialScale)
  model_x <- RMgauss(var=SD_x^2, scale=SpatialScale)
  model_eta <- RMgauss(var=SD_eta^2, scale=SpatialScale)
  
  # Realization from GRF
  # delta_s -- spatial variation that affects both density and sampling intensity
  # x_s -- spatial variation in 
  delta_s = RFsimulate(model=model_delta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
  x_s = RFsimulate(model=model_x, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
  eta_s = RFsimulate(model=model_eta, x=loc_s[,'x'], y=loc_s[,'y'])@data[,1]
  
  # Total abundance
  Ztrue_s = exp( beta0 + betax*x_s + delta_s )
  Ztrue_s = rpois(n_cells,Ztrue_s)     #added to v5, 11/9/2015

  Data[1,irow:(irow+n_cells-1)]=(Ztrue_s-mean(Ztrue_s))/mean(Ztrue_s)  
  #simulate sampling
  for(SimI in 1:length(b_set)){
      # Samping intensity
      b=b_set[SimI]
      R_s = exp( betax_prob*x_s + eta_s + b*delta_s )
      R_s = R_s / sum(R_s)
      Data[SimI+1,irow:(irow+n_cells-1)]=R_s
  }
  irow=irow+n_cells
}
  
library(mgcv)
gam_data=data.frame(t(Data))
colnames(gam_data)=c("Resid","R0","R1","R5")
b0=gam(R0~s(Resid),data=gam_data)
b1=gam(R1~s(Resid),data=gam_data)
b5=gam(R5~s(Resid),data=gam_data)

#make predictions over same range of abundance residuals
New_data=data.frame(Resid=-1+c(1:100)*(6/100))
b0.pred=predict(b0,newdata=New_data)
b1.pred=predict(b1,newdata=New_data)
b5.pred=predict(b5,newdata=New_data)


Plot.data=data.frame(N.residual=rep(unlist(New_data),3),Prob.sampled=c(b0.pred,b1.pred,b5.pred),b=rep(factor(c(0,1,5)),each=nrow(New_data)))
Density=density(gam_data$Resid,from=-1,to=5,n=100)
Plot.data$KDE=rep(Density$y,3)*max(Plot.data$Prob.sampled)
library(ggplot2)
my_plot=ggplot(data=Plot.data,aes(x=N.residual))+geom_line(aes(y=Prob.sampled,linetype=b,col=b),size=1.3)+geom_line(aes(y=KDE),size=1.3)
my_plot=my_plot + theme(text=element_text(size=20))+theme(axis.text.y=element_text(size=14))+labs(x = "Standardized abundance residual", y=expression(paste("Sampling probability (",R["i"],")")))

my_plot

pdf("sim_R.pdf")
my_plot
dev.off()
  