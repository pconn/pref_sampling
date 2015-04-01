set.seed(123)
library(sp)
library(gstat)
library(fields)
library(PCspatreg)
xy <- expand.grid(1:100, 1:100)
names(xy) <- c('x','y')
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vgm(psill=1, range=20, nugget=0.01, model='Gau'), nmax=20)
simData <- predict(g.dummy, newdata=xy, nsim=4)
simData$sim4 = simData$sim4-mean(simData$sim4)

simData$resp = 3 + 2*simData$sim1 + 0.1*simData$sim3 + (simData$sim4 + rnorm(1000,0,1))

idx = sort(sample(1:10000, 100, replace=FALSE))
sampData = simData[idx,]

knots = expand.grid(x=seq(-10, 110, 15), y=seq(-10, 110, 15))

plot(knots$x, knots$y, asp=1, type="n")
image(x=1:100, y=1:100, z=t(matrix(simData$resp,100,100)[100:1,]), asp=TRUE, add=TRUE)
points(knots$x, knots$y)
points(sampData$x, sampData$y, pch=16)

K = zapsmall(exp(-0.5*rdist(simData[,c("x","y")], knots)^2/15^2)/(15*sqrt(2*pi)))
K = sweep(K, 2, colSums(K),"/")/81
X = as.matrix(cbind(1, simData[,3:5]))

betaMean = c(mean(simData$resp[idx]), rep(0,3))
betaSig=10*solve(crossprod(X[idx,]))


postSamp = PCspatregMCMC(y=sampData$resp, X=X, K=K, sampleIndex=idx, betaMean=betaMean, betaSig=betaSig, 
                         phiScale=10, burn=1000, iter=10000)

library(coda)
plot(mcmc(postSamp$beta))
plot(mcmc(abs(postSamp$phi)))
plot(mcmc(postSamp$sigma))



simData$pred = colMeans(postSamp$pred)
plot(knots$x, knots$y, asp=1, type="n")
image(x=1:100, y=1:100, z=t(matrix(simData$pred,100,100)[100:1,]), asp=TRUE, add=TRUE)
points(knots$x, knots$y)
points(sampData$x, sampData$y, pch=16)

