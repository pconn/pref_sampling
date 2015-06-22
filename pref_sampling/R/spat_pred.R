#' function to perform Bayesian analysis of log Gaussian Cox processes and to make spatial predictions 
#' @param formula A formula object specifying the model for the fixed effects
#' @param Data   A SpatialPolygonsDataFrame holding covariate values for all cells one wants to make predictions on (including sampled cells)
#' @param Effort A list including "Counts" which gives the response variable for surveyed cells, and "Mapping" which is a vector indicating which sample unit is associated with which count
#' @param Area.adjust   A vector allowing for differences in suitable habitat for each cell.  Can be used for different grid cell sizes or different
#'        proportions of suitable habitat (e.g., 1.0 = 100 percent of habitat is suitable for the focal species)
#' @param Offset	A vector giving the fraction of each grid cell in Mapping that was surveyed
#' @param spat.mod Controls whether (and what type) of spatial autocorrelation is implemented (spatmod = 0: no spatial autocorrelation; spatmod=1: process convolution; spatmod=2: srr ICAR) 
#'        Note that if spatmod=1, the argument Knots needs to be specified; if spatmod=2, Adj must be specified
#' @param Assoc   Association matrix for all grid cells in Data (needed for ICAR models only)
#' @param K   A (# cells)x(# knots) matrix - the 'K' matrix for process convolution modeling.  
#' @param Names.gam A character vector giving the column names of Data that will be modeled as smooth effects
#' @param Knots.gam A list vector (length is # of smooth params), where each list element is a vector providing locations of knots (one location for each knot)
#' @param Control A list giving MCMC controls and additional model options
#'  "iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'  "srr.tol": Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation (default is 0.5)
#'  "predict" If TRUE (default), calculate posterior predictions across the entire grid
#'  "MH.nu" A vector providing continuous Uniform half-range values for Metropolis-Hastings proposals
#'  "adapt" If true, adapts tuning parameters for nu updates; e.g., one could run a chain with adapt=TRUE, and use adapted MH tuning parameters in a longer analysis (default FALSE)
#'  "fix.tau.epsilon" If TRUE, fixes tau.epsilon to 100
#'  "Kern.gam.se" A vector setting gam se for kernels (length is number of smooth params)
#' @param Prior.pars  An optional list giving prior parameters; these include
#'  "beta.tau" precision for Gaussian prior on regression parameters (default 0.1)
#'  "a.eps" alpha parameter for tau_epsilon~Gamma(alpha,beta) (default = 1.0)
#'  "b.eps" beta parameter for tau_epsilon~Gamma(alpha,beta) (default = 0.01)
#'   If provided, all values for Prior.pars need to be specified
#' @param Precision.pars If provided, this list gives tau.epsilon and tau.eta, and these parameters are both fixed during estimation
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#'  gIVH: A binary vector, with length equal to the number of survey units, describing whether predictions for each survey unit are (=1) or are not (=0) within the gIVH
#'  K.gam: A ``design matrix" associated with smooth effects, should the GAM procedure be employed
#' @export
#' @import Matrix
#' @keywords areal model, data augmentation, mcmc, spatial prediction
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
#' @examples print("Later!")
spat_pred<-function(formula,Data,Effort,spat.mod=0,Offset,Area.adjust,Control,Assoc=NULL,K=NULL,Names.gam=NULL,Knots.gam=NULL,Prior.pars=NULL,Precision.pars=NULL){
  #do a little checking
  if(spat.mod==1 & is.null(K)==TRUE)cat('ERROR: Knots Spatial Points Object (see sp library) must be supplied for spatmod=1')  
  if(spat.mod==2 & is.null(Assoc)==TRUE)cat('ERROR: Association matrix must be supplied for spatmod=2')

  if(is.null(Prior.pars)){
    Prior.pars=list(beta.tau=0.01,
    a.eps=1,
    b.eps=0.01,
    a.eta=1,
    b.eta=0.01,
    a.gam=1,
    b.gam=0.01)
  }
  Count=Effort$Counts
  Mapping=Effort$Mapping
  Offset=Offset*Area.adjust[Mapping]
  Log.offset=log(Offset)
  Log.area.adjust=log(Area.adjust)
  S=nrow(Data@data) #number of cells
  n=length(Mapping) #number of cells with data
  mcmc.length=(Control$iter-Control$burnin)/Control$thin
  
  #setup design matrix

  X.pred=model.matrix(formula,data=Data@data)
  X.obs=model.matrix(formula,data=Data@data[Mapping,])

  n.beta=ncol(X.obs)
  
  #check to make sure DMs same size
  if(ncol(X.pred)!=n.beta)cat('ERROR: prediction and observation design matrices have different number of columns. Check to make sure all levels of categorical variables are represented in sampled cells')

  Which.not.sampled=c(1:S)
  Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]
  n.no=S-n
  X.no=as.matrix(X.pred[Which.not.sampled,],ncol=ncol(X.obs))

  XpXinv=solve(crossprod(X.obs)) 
  XpXinvXp=XpXinv%*%t(X.obs)


  #GAM stuff
  if(is.null(Names.gam)==FALSE){
    n.smooths=length(Names.gam)
    Tmp=vector("list",n.smooths)
    #formulate K matrix for each smooth covariate
    for(icov in 1:n.smooths){
      Tmp[[icov]]=matrix(0,S,length(Knots.gam[[icov]])) 
      for(iknot in 1:length(Knots.gam[[icov]])){
        Tmp[[icov]][,iknot]=dnorm(Data@data[,Names.gam[icov]],Knots.gam[[icov]][iknot],Control$Kern.gam.sd[icov])
      }
      #Tmp[[icov]]=Tmp[[icov]]/rowSums(Tmp[[icov]])  #this produces singular KpK if >1 smooth term
    }
    K.gam=Tmp[[1]]
    if(n.smooths>1){
      for(icov in 2:n.smooths)K.gam=cbind(K.gam,Tmp[[icov]])
    }
    #K.gam=K.gam/apply(K.gam,1,'sum')  don't do this - multimodal GAMs!
    n.alpha=ncol(K.gam)
    Alpha=rep(0,n.alpha)
    cross.K.gam=crossprod(K.gam[Mapping,])
    Alpha.mc=matrix(0,n.alpha,mcmc.length)
    #KpKinv.gam=solve(cross.K.gam)
  }
  Gamma=rep(0,S)  #combined smooth effects by site
  
  #provide initial estimate of Nu
  Nu=rep(0,S)
  Nu[Mapping]=Count/Offset
  Nu[which(Nu==0)]=min(Nu[which(Nu>0)])
  #Nu[Which.not.sampled]=mean(Nu[Mapping])
  Nu=log(Nu)+rnorm(S,0,0.05)
  Beta=rep(0,n.beta)
  Beta[1]=mean(Nu) #set intercept to reasonable starting value
  Beta.mc=matrix(0,n.beta,mcmc.length)
  tau.epsilon=100
  tau.eta=100
  if(is.null(Precision.pars)==FALSE){
    tau.epsilon=Precision.pars$tau.epsilon
    tau.eta=Precision.pars$tau.eta
  }
  tau.alpha=100
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.eta.mc=Tau.epsilon.mc
  Tau.alpha.mc=Tau.epsilon.mc #for gam
  Pred.mc=matrix(0,S,mcmc.length)
  Eta=rep(0,S)
  Accept=rep(0,n)
  Accept.old=Accept
  
  if(spat.mod==1){ #define some matrices, etc. needed for PC model
    n.knots=ncol(K)
    K.t=t(K)
    cross.K=crossprod(K)
    Theta=rep(0,n.knots)
    KpKinv=solve(cross.K) # add a little noise to make matrix non-singular
    KpKinvKp=XpXinv%*%t(X.obs)
    n.theta=n.knots
  }
  
  if(spat.mod==2){ #define some matrices, etc. needed for spatial updates under srr model
    Q=-Assoc
    diag(Q)=apply(Assoc,2,'sum')
    Q=Matrix(Q)	
    L.t=XpXinv
    L=L.t
    Qt=L.t
    cross.L=L.t
    Theta=L.t
    P.c=diag(S)-X.pred%*%solve(crossprod(X.pred),t(X.pred))
    Omega=Matrix((P.c%*%Assoc%*%P.c)*(S/sum(Assoc)))
    Eigen=eigen(Omega)
    if(Control$srr.tol>1){ #in this case, the number of eigenvalues is preselected
      L.t=Eigen$vectors[,1:Control$srr.tol]      
    }
    else{ #in this case, compute the number of eigenvalues > the threshold
      if(max(Eigen$values)<Control$srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
      Ind=which(Eigen$values>Control$srr.tol)
      L.t=Eigen$vectors[,Ind]
      cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
    }
    L=t(L.t)
    Qt=L%*%Q%*%L.t
    cross.L=L%*%L.t
    n.theta=nrow(Qt)
    Theta=rnorm(n.theta,0,sqrt(1/tau.eta))
  }
  
  if(spat.mod>0)Theta.mc=matrix(0,n.theta,mcmc.length)
     
  if(is.null(Names.gam)==FALSE & n.beta==1){
    Beta[1]=0
    n.beta=0
  }
    
  set.seed(12345)
  
  for(iiter in 1:Control$iter){
    if(iiter%%1000==0)cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))
    #update nu (sampled cells)
    Mu=X.obs%*%Beta+Eta[Mapping]+Gamma[Mapping]
    sd=sqrt(1/tau.epsilon)
    full.cond.old=dnorm(Nu[Mapping],Mu,sd,log=1)+dpois(Count,exp(Log.offset+Nu[Mapping]),log=1)
    Prop=Nu[Mapping]+runif(n,-Control$MH.nu,Control$MH.nu)
    full.cond.new=dnorm(Prop,Mu,sd,log=1)+dpois(Count,exp(Log.offset+Prop),log=1)
    I.accept=(runif(n)<exp(full.cond.new-full.cond.old))
    Nu[Mapping][I.accept==1]=Prop[I.accept==1]
    Accept=Accept+I.accept
    
    #simulate nu for unsampled cells at each iteration IF using CAR model (otherwise just do for iterations where MCMC is stored)
    if(spat.mod==2)Nu[Which.not.sampled]=rnorm(n.no,X.no%*%Beta+Eta[Which.not.sampled]+Gamma[Which.not.sampled],sqrt(1/tau.epsilon))
    #Nu=Nu.true
    #Nu=log(Data@data[["Lambda"]])

    #update beta
    if(n.beta>0)Beta=t(rmvnorm(1,XpXinvXp%*%(Nu[Mapping]-Eta[Mapping]-Gamma[Mapping]),XpXinv/(tau.epsilon+Prior.pars$beta.tau)))
    
    #update precision for exchangeable errors
    if(Control$fix.tau.epsilon==FALSE & is.null(Precision.pars)==TRUE){
      Mu=X.obs%*%Beta+Eta[Mapping]+Gamma[Mapping]
      Diff=Nu[Mapping]-Mu
      tau.epsilon <- rgamma(1,n/2 + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
    }
    
    #update SRR parameters, spatial random effects if applicable
    if(spat.mod==2){
      Dat.minus.Exp=Nu-X.pred%*%Beta-Gamma
      V.eta.inv <- cross.L*tau.epsilon + tau.eta*Qt
      M.eta <- solve(V.eta.inv, tau.epsilon*L%*%Dat.minus.Exp)
      Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.theta,0,1))
      Eta=as.numeric(L.t%*%Theta)					
      #update tau.eta
      if(is.null(Precision.pars)==TRUE)tau.eta <- rgamma(1, n.theta*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta, Qt %*% Theta)*0.5) + Prior.pars$b.eta)    
    }
    
    if(spat.mod==1){ #update process convolution parameters/effects
      Dat.minus.Exp=Nu-X.pred%*%Beta-Gamma
      V.eta.inv <- cross.K*tau.epsilon + diag(tau.eta,nrow=n.knots)
      M.eta<- solve(V.eta.inv, tau.epsilon*K.t%*%Dat.minus.Exp)
      Theta <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.knots,0,1))
      Eta=as.numeric(K%*%Theta)
      #update tau.eta
      if(is.null(Precision.pars)==TRUE)tau.eta <- rgamma(1, n.knots*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Theta))*0.5 + Prior.pars$b.eta)          
    }
    
    #update GAM parameters
    if(is.null(Names.gam)==FALSE){
      Dat.minus.Exp=Nu[Mapping]-X.obs%*%Beta-Eta[Mapping]
      #Alpha=rmvnorm(1,KpKinv.gam%*%crossprod(K.gam[Mapping,],Dat.minus.Exp),KpKinv.gam/(tau.epsilon+tau.alpha))
      V.alpha.inv <- cross.K.gam*tau.epsilon + diag(tau.alpha,nrow=n.alpha)
      #V.alpha.inv <- cross.K.gam*(tau.epsilon +tau.alpha)
      M.alpha <- solve(V.alpha.inv,tau.epsilon*t(K.gam[Mapping,])%*%Dat.minus.Exp)
      Alpha <- M.alpha + solve(chol(as.matrix(V.alpha.inv)), rnorm(n.alpha,0,1))
      Gamma=as.numeric(K.gam%*%Alpha)
      #Gamma=as.numeric(K.gam%*%t(Alpha))
      #update tau.alpha
      tau.alpha <- rgamma(1,n.alpha*0.5 + Prior.pars$a.gam,as.numeric(0.5*crossprod(Alpha)) + Prior.pars$b.gam)    
    }
    
    #adapt proposal distributions if Control$adapt=TRUE
    if(Control$adapt==TRUE & iiter%%100==0){
      Diff=Accept-Accept.old
      for(inu in 1:n){
        if(Diff[inu]<30)Control$MH.nu[inu]=Control$MH.nu[inu]*.95
        if(Diff[inu]>40)Control$MH.nu[inu]=Control$MH.nu[inu]*1.053
      }
      Accept.old=Accept
    }
    
    #store MCMC values when appropriate (including predictions)
    if(iiter>Control$burnin & iiter%%Control$thin==0){
      Beta.mc[,(iiter-Control$burnin)/Control$thin]=Beta
      Tau.epsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.epsilon
      if(is.null(Names.gam)==FALSE){
        Tau.alpha.mc[(iiter-Control$burnin)/Control$thin]=tau.alpha
        Alpha.mc[,(iiter-Control$burnin)/Control$thin]=Alpha
      }
      if(spat.mod>0){
        Tau.eta.mc[(iiter-Control$burnin)/Control$thin]=tau.eta
        Theta.mc[,(iiter-Control$burnin)/Control$thin]=as.numeric(Theta)
      }
      if(Control$predict==TRUE){ #make predictions
        #simulate nu if not an ICAR model
        if(spat.mod<2)Nu[Which.not.sampled]=rnorm(n.no,X.no%*%Beta+Gamma[Which.not.sampled]+Eta[Which.not.sampled],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        Pred.mc[,(iiter-Control$burnin)/Control$thin]=rpois(S,exp(Log.area.adjust+Nu))    
      }
    }
  } 
  Beta.var=t(cov(t(Beta.mc)))
  Mu.var=X.pred%*%tcrossprod(Beta.var,X.pred)
  Mu.mat=Pred.mc
  for(i in 1:mcmc.length)Mu.mat[,i]=X.pred%*%Beta.mc[,i]
  
  Out=list(MCMC=list(Beta=Beta.mc,tau.epsilon=Tau.epsilon.mc,Pred=Pred.mc),Accept=Accept,Control=Control)
  if(spat.mod>=1){
    Out$MCMC$tau.eta=Tau.eta.mc
    Out$Eta=Eta
    Out$Sigma.beta=array(0,dim=c(dim(XpXinv),length(Tau.epsilon.mc)))
    for(i in 1:length(Tau.epsilon.mc))Out$Sigma.beta[,,i]=XpXinv*(Tau.epsilon.mc[i]+Prior.pars$beta.tau)^-1
  }
  if(spat.mod==1){
    Out$K.sp=K
    Out$Sigma.alpha=array(0,dim=c(dim(cross.K),length(Tau.epsilon.mc)))
    for(i in 1:length(Tau.epsilon.mc))Out$Sigma.alpha[,,i]=solve(cross.K*Tau.epsilon.mc[i] + diag(Tau.eta.mc[i],nrow=n.knots))                    
    Mu.var=Mu.var+K%*%tcrossprod(t(cov(t(Theta.mc))),K)
    for(i in 1:mcmc.length)Mu.mat[,i]=Mu.mat[,i]+K%*%Theta.mc[,i]
  }
  if(spat.mod==2){
    Out$K.sp=L.t
    Out$Sigma.alpha=array(0,dim=c(dim(cross.L),length(Tau.epsilon.mc)))
    for(i in 1:length(Tau.epsilon.mc))Out$Sigma.alpha[,,i]=as.matrix(solve(cross.L*Tau.epsilon.mc[i]+Tau.eta.mc[i]*Qt))
    Mu.var=Mu.var+L.t%*%t(cov(t(Theta.mc)))%*%L
    for(i in 1:mcmc.length)Mu.mat[,i]=Mu.mat[,i]+L.t%*%Theta.mc[,i]
  }
  if(is.null(Names.gam)==FALSE){
    Out$MCMC$tau.alpha=Tau.alpha.mc
    Out$MCMC$Alpha=Alpha.mc
    Out$K.gam=K.gam
    Mu.var=Mu.var+K.gam%*%tcrossprod(t(cov(t(Alpha.mc))),K.gam)
    for(i in 1:mcmc.length)Mu.mat[,i]=Mu.mat[,i]+K.gam%*%Alpha.mc[,i]
  }
  Mu=apply(Mu.mat,1,'median')
  Lambda.var=(exp(Mu))^2*diag(Mu.var)
  max.obs=max(Lambda.var[Mapping])
  gIVH=rep(1,S)
  which.gt.max=which(Lambda.var>max.obs)
  if(length(which.gt.max)>0)gIVH[which.gt.max]=0
  max.obs=max(diag(Mu.var)[Mapping])
  Out$gIVH=gIVH
  Out
}

