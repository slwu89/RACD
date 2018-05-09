library(RACD)
library(RACDaux)
library(tidyverse)
library(ggplot2)
library(ggforce)

theta = RACDaux::RACD_Parameters()
init = RACDaux::RACD_Setup(theta)
## Not run: 
library(tidyverse)
outfile = "/Users/pandagoneamok/Desktop/tmp.csv"
RACD_Simulation(365,theta,init$humans,init$houses,123,outfile)
state = RACDaux::RACD_StateVector(outfile)
state %>% as.tibble %>% gather(state,value,-time) %>% ggplot(aes(x=time,y=value,color=state)) + geom_line() + theme_bw()

calcDistance <- function(x, y, breedingSite) {
  # x : x-coordinate
  # y : y-coordinate
  # breedingSite : list of breeding site parameters (x, y, sigma)
  sqrt(sum((x-breedingSite$x)^2, (y-breedingSite$y)^2))
}

initRiskGrid <- function(breedingSiteList, n) {
  # breedingSiteList : list of breeding sites (x, y, sigma)
  # n : number of points to split x- and y-axes into
  riskGrid <- expand.grid((1:n)/n, (1:n)/n)
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

plotInitCond <- function(initCond, n) {
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
    lapply(houseState,
           function(house) {
             geom_arc_bar(data=house, 
                          aes(x0=x, y0=y, r0=0, r=pieR, amount=Freq, colour=state), 
                          stat="pie")
           })
  return(riskPlot)
}

plotInitCond(init, n=100)

## plot immunity functions
# get fixed parameters (supplemental)

## table of model parameters
# ?RACD_Parameters

## plot simPoints.R