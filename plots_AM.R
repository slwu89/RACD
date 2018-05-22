library(RACD)
library(RACDaux)
library(tidyverse)
library(ggforce)
library(gridExtra)

## risk contour map with breeding sites and houses

calcDistance <- function(x, y, breedingSite) {
  # x : x-coordinate
  # y : y-coordinate
  # breedingSite : list of breeding site parameters (x, y, sigma)
  sqrt(sum((x-breedingSite$x)^2, (y-breedingSite$y)^2))
}

initRiskGrid <- function(breedingSiteList, n) {
  # breedingSiteList : list of breeding sites (x, y, sigma)
  # n : number of points to split x- and y-axes into
  axes <- seq(from=-0.05, to=1.05, length.out=n)
  riskGrid <- expand.grid(axes, axes)
  colnames(riskGrid) <- c("x", "y")
  riskGrid$risk <- sapply(1:nrow(riskGrid),
                          function(i) {
                            sum(sapply(breedingSiteList,
                                       function(site) {
                                         dnorm(calcDistance(riskGrid[i,1], riskGrid[i,2], site),
                                               sd=site$sigma)
                                       }))
                          })
  return(riskGrid)
}

plotInitCond <- function(initCond, n, pie=T) {
  # initCond : list of humans, breedingSites, and houses
  # humans : list of human parameters
  # breedingSites : list of breedingSite parameters (x, y, sigma)
  # houses : list of house parameters (x, y, size, psi)
  
  # compute risk grid
  riskGrid <- initRiskGrid(initCond$breedingSites, n)
  
  # unlist breedingSites, houses, and humans
  breedingSites <- data.frame(t(sapply(initCond$breedingSites, unlist)))
  houses <- data.frame(t(sapply(initCond$houses, unlist)))
  humans <- data.frame(t(sapply(initCond$humans, unlist)))
  
  # compute proportion of states at each house
  houseState <- data.frame(table(humans$house, humans$state))
  colnames(houseState)[1:2] <- c("house", "state")
  houseState <- cbind(houseState, houses[match(houseState$house, rownames(houses)),1:2])
  pieR <- diff(range(houseState$x))/50
  houseState <- split(houseState, houseState$house)
  
  # plot
  riskPlot <- ggplot() +
    geom_tile(data=riskGrid, aes(x=x, y=y, fill=risk)) +
    stat_contour(data=riskGrid, aes(x=x, y=y, z=risk)) +
    geom_point(data=breedingSites, aes(x=x, y=y), color="darkred") +
    # geom_point(data=houses, aes(x=x, y=y), color="chartreuse3") +
    theme_bw() + xlab("") + ylab("")
  if(pie) {
    riskPlot <- riskPlot +
      lapply(houseState,
             function(house) {
               geom_arc_bar(data=house, 
                            aes(x0=x, y0=y, r0=0, r=pieR, amount=Freq, colour=state), 
                            stat="pie")
             })
  }
  return(riskPlot)
}

theta = RACDaux::RACD_Parameters()
init = RACDaux::RACD_Setup(theta)
plotInitCond(init, n=100)
plotInitCond(init, n=100, pie=F)

grid.arrange(plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             ncol=2)

## plot simPoints.R
# xy_pts : tibble with point type and (x,y) coordinates

plot_xy_pts <- function(xy_pts, house=T, breedingSite=T) {
  # xy_pts : tibble with point type (habitat/house) and (x,y) coordinates
  # house : logical, whether to plot houses
  # breedingSites : logical, whether to plot mosquito breeding sites
  
  # split xy_pts into houses and habitats
  xy_h <- subset(xy_pts, type=="houses")
  xy_b <- subset(xy_pts, type=="habitats")
  
  # compute sigma for each habitat
  xy_b$sigma <- sapply(1:nrow(xy_b),
                       function(i) {
                         min(as.matrix(dist(x=rbind(xy_b[i,2:3],xy_h[,2:3])))[1,2:(nrow(xy_h)+1)])
                       })
  
  # compute risk grid
  axes <- seq(from=-0.1, to=2.1, length.out=100)
  riskGrid <- expand.grid(axes, axes)
  colnames(riskGrid) <- c("x", "y")
  riskGrid$risk <- sapply(1:nrow(riskGrid),
                          function(i) {
                            risks <- sapply(1:nrow(xy_b),
                                              function(j) {
                                                dist <- sqrt(sum((riskGrid$x[i]-xy_b$x[j])^2,
                                                                 (riskGrid$y[i]-xy_b$y[j])^2))
                                                dnorm(dist, sd=xy_b$sigma[j])
                                              })
                            sum(risks)
                          })
  
  # plot
  riskPlot <- ggplot() +
    geom_tile(data=riskGrid, aes(x=x, y=y, fill=risk)) +
    stat_contour(data=riskGrid, aes(x=x, y=y, z=risk)) +
    # scale_fill_gradient2(low="black", mid="dodgerblue4", high="lightsteelblue2") +
    theme_bw() + xlab("") + ylab("")
  if(house) {
    riskPlot <- riskPlot + geom_point(data=xy_h, aes(x=x, y=y), color="chartreuse3")
  }
  if(breedingSite) {
    riskPlot <- riskPlot + geom_point(data=xy_b, aes(x=x, y=y), color="darkred", 
                                      alpha=0.5, size=0.5)
  }
  
  return(riskPlot)
}

plot_xy_pts(xy_pts)

## plot immunity functions
# get fixed parameters (supplemental)
