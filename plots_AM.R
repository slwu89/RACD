library(RACD)
library(RACDaux)
library(tidyverse)
library(ggforce)
library(gridExtra)

theta = RACDaux::RACD_Parameters()
init = RACDaux::RACD_Setup(theta)

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
  riskGrid <- initRiskGrid(initCond$breedingSites, n)
  breedingSites <- data.frame(t(sapply(initCond$breedingSites, unlist)))
  houses <- data.frame(t(sapply(initCond$houses, unlist)))
  humans <- data.frame(t(sapply(initCond$humans, unlist)))
  houseState <- data.frame(table(humans$house, humans$state))
  colnames(houseState)[1:2] <- c("house", "state")
  houseState <- cbind(houseState, houses[match(houseState$house, rownames(houses)),1:2])
  pieR <- diff(range(houseState$x))/50
  houseState <- split(houseState, houseState$house)
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

plotInitCond(init, n=100)
plotInitCond(init, n=100, pie=F)

grid.arrange(plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             plotInitCond(RACDaux::RACD_Setup(theta), n=100),
             ncol=2)

## plot immunity functions
# get fixed parameters (supplemental)

## table of model parameters
# ?RACD_Parameters

## plot simPoints.R