library(spatstat)
library(mvtnorm)
library(RACD)
library(RACDaux)

################################################################################
# simulate pointset
###############################################################################

# spatially varying intensity function for houses (houses cluster together like gaussians)
win <- square(1)
sigma1 <- diag(2)*0.005
sigma2 <- diag(2)*0.005
mean1 <- c(0.75,0.75)
mean2 <- c(0.25,0.25)
w = 2/3
fun_house <- function(x,y,mean1,mean2,sigma1,sigma2,w){
  (w*dmvnorm(x =c(x,y),mean =mean1,sigma = sigma1,log = FALSE)) + 
  ((1-w)*dmvnorm(x = c(x,y),mean = mean2,sigma = sigma2,log = FALSE))
}
fun_house <- Vectorize(FUN = fun_house,vectorize.args = c("x","y"))
mu_house <- as.im(fun_house,W = win,mean1=mean1,mean2=mean2,sigma1=sigma1,sigma2=sigma2,w=w)

# spatially varying intensity function for breeding sites (they are more dispersed, like a t-distribution)


# we can choose to simluate how houses arise as a doubly stochastic process
# in this case the underlying intensity is a random field, here modeled as
# log(intensity) ~ GP(.)
rf_intensity <- rLGCP(model = "gauss",mu = log(mu_house),var = 0.15,scale = 0.1,win = win)
rf_intensity <- attr(rf_intensity, "Lambda")
par(mfrow=c(1,2))
plot(rf_intensity)
plot(mu_house)
par(mfrow=c(1,1))

# simulate houses
n_house <- 100

# use Metropolis-Hastings simulation for point process because we need to condition on number
# of houses in the inhomogeneous  point process
# pois_mod <- list(cif="poisson",par=list(beta=1),w=win,trend=rf_intensity)
hc_mod <- list(cif="hardcore",par=list(beta=1,hc=0.025),w=win,trend=rf_intensity)
rmh_control <- list(p=1)
rmh_start <- list(n.start=n_house)
xy_h <- rmh(model = hc_mod,control = rmh_control,start = rmh_start)
plot(xy_h)

# simulate breeding sites, assume breeding sites cluster around houses
