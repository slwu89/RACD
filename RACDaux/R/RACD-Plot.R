###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   Plot Landscape
#
###############################################################################

# if (i%%1==0) { # Visualize cases daily
#   # 1. Plot risk surface on a 1x1 grid:
#   # First, calculate risk surface for the map grid:
#   riskMap <- matrix(rep(0,101*101),101)
#   for (j in 1:101) {
#     longMap <- (j-1)/100 # Convert index to 0-1 map scale
#     for (k in 1:101) {
#       latMap <- (k-1)/100 # Convert index to 0-1 map scale
#       for (l in 1:numBreedingSites) {
#         # Compute cumulative risk due to all of the breeding sites:
#         riskMap[j,k] <- ( riskMap[j,k] + (1/(sigmaBreedingSite[l]*sqrt(2*pi)))
#                 * exp(-((longMap-longBreedingSite[l])^2 +
#                                      (latMap-latBreedingSite[l])^2)
#                                     / (2*sigmaBreedingSite[l]^2)))
#       }
#     }
#   }
#   numContours <- 9
#   contour(x = seq(0, 1, length.out = nrow(riskMap)),
#       y = seq(0, 1, length.out = ncol(riskMap)),
#       z= riskMap, drawlabels=FALSE, nlevels=numContours, lwd=2,
#       col=brewer.pal(numContours, "Blues"), axes=FALSE)
#
#   # 2. Plot breeding sites
#   # for (j in 1:numBreedingSites) {
#   #	points(longBreedingSite[j], latBreedingSite[j], col="red", lwd=3)
#   # }
#
#   # 3. Plot houses
#   houseEdgeSize <- 0.015
#   for (j in 1:numHouses) {
#     polygon(x=c(longHouse[j]+houseEdgeSize/2, longHouse[j]+houseEdgeSize/2, longHouse[j]-houseEdgeSize/2, longHouse[j]-houseEdgeSize/2, longHouse[j]),
#             y=c(latHouse[j]+houseEdgeSize/2, latHouse[j]-houseEdgeSize/2, latHouse[j]-houseEdgeSize/2, latHouse[j]+houseEdgeSize/2, latHouse[j]+houseEdgeSize),
#         col="brown", border="black", lwd=2)
#   }
#
#   # 3. Plot T (green), D (red), A (orange), U (yellow) as they come
#   #    up (around house):
#   radCases <- 0.025
#   for (j in 1:numHouses) {
#     isTDAUHouseJ <- which( ( sapply(indiv, function(x) x$house) == j) &
#                ( ( sapply(indiv, function(x) x$state) == "T") |
#              ( sapply(indiv, function(x) x$state) == "D") |
#              ( sapply(indiv, function(x) x$state) == "A") |
#              ( sapply(indiv, function(x) x$state) == "U") ) )
#     numCasesHouseJ <- length(isTDAUHouseJ)
#     if (numCasesHouseJ > 0) {
#       for (k in 1:numCasesHouseJ) {
#         if (indiv[[isTDAUHouseJ[k]]]$state == "T") {
#           points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="green", lwd=2, pch=16)
#         } else if (indiv[[isTDAUHouseJ[k]]]$state == "D") {
#           points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="red", lwd=2, pch=16)
#         } else if (indiv[[isTDAUHouseJ[k]]]$state == "A") {
#           points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="orange", lwd=2, pch=16)
#         } else if (indiv[[isTDAUHouseJ[k]]]$state == "U") {
#           points(longHouse[j] + radCases*sin((k-1)*2*pi/numCasesHouseJ), latHouse[j] - radCases*cos((k-1)*2*pi/numCasesHouseJ), col="yellow", lwd=2, pch=16)
#         }
#       }
#     }
#   }
# }
