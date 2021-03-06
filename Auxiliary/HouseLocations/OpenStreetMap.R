###############################################################################
###############################################################################
# imagine ascii art
# that will happen
# when I
# have time.
###############################################################################
###############################################################################
# safety
invisible(rm(list=ls()));invisible(gc());

# load libraries
library(osmar)
library(sp)

#######################################
# functions from package that I needed more options in.
# these cover the OSMAR functions in the global environment
merge_ways_nodes <- function(ways, nodes) {
  colnames(ways) <- sprintf("w%s", colnames(ways))
  colnames(nodes) <- sprintf("n%s", colnames(nodes))
  
  m <- match(ways$wref, nodes$nid)
  
  dat <- cbind(ways, nodes[m, ])
  # dat <- na.omit(dat)
  dat <- dat[!is.na(dat$nlat), ]
  
  dat$nid <- NULL
  colnames(dat) <- substring(colnames(dat), 2)
  
  dat
}

plot_ways <- function(x, add = FALSE, xlab = "lon", ylab = "lat",main = "Location",
                      ...) {
  #stopifnot(is_osmar(x))
  
  dat <- merge_ways_nodes(x$ways[[3]], x$nodes[[1]])
  
  rlat <- range(dat$lat, na.rm = TRUE)
  rlon <- range(dat$lon, na.rm = TRUE)
  
  dat <- split(dat, dat$id)
  
  if ( !add ) {
    plot(1, type = "n", xlim = rlon, ylim = rlat,
         xlab = xlab, ylab = ylab, main = main)
  }
  
  for ( coord in dat ) {
    #if ( nrow(coord) >= 2 ) {
    lines(coord[, c("lon", "lat")], ...)
    #}
  }
}
###############################################################################
###############################################################################
###############################################################################
# 1a) Download osm file for desired country/region/area
#    https://download.geofabrik.de/
# 1b) Or, use open street map to pick an area and download. This is my preferred 
#     https://www.openstreetmap.org

# THESE ROUTINES CAN'T HANDLE LARGE NUMBERS OF NODES!!! NO FILES OVER ~2.3 GB, AND 
#  IT CHOKES ON THINGS AS LARGE AS 750MB. NEED TO STICK WITH SMALL CITIES FOR NOW.

#######################################

# Import data, This can take awhile
OSMData <- get_osm(x = complete_file(), source = osmsource_file(file = "./stJohnsUniversityTanzania.osm"))

title <- "st John's University Tanzania"

###############################################################################
###############################################################################
# subset roads
road_ids <- OSMData$ways$tags$id[ OSMData$ways$tags$k == "highway" ]
road_ids <- find_down(object = OSMData, ids = way(road_ids))
roads <- subset(x = OSMData, ids = road_ids)

# subset buildings
bg_ids <- OSMData$way$tags$id[ OSMData$ways$tags$k == "building" ]
bg_ids <- find_down(OSMData, way(bg_ids))
bg <- subset(OSMData, ids = bg_ids)

# convert to sp object, calculate centroids. This is the only reason for the sp package,
#  and the conversion is slow as fuck. Need to write a better centroids calculator.
bg_poly <- as_sp(bg, "polygons")
bg_centroids <- sp::coordinates(bg_poly)

###############################################################################
###############################################################################
# scale and setup output
#horiz <- (left - right) / (top - bottom)

horiz <- diff(range(bg$nodes$attrs$lon)) / diff(range(bg$nodes$attrs$lat))
png(paste0("./", gsub(pattern = " ", replacement = "_", x = title, fixed = TRUE), ".png"),
    width=2304*horiz, height=2304, units="px")

# plot
par(las = 1, font.lab = 2, font.axis = 2, font.main = 2, cex.main = 4,
    cex.lab = 3, cex.axis = 2, mar=c(7,8,5,2), mgp = c(5,1,0))

plot_ways(x = bg,  main = title, ylab = "Latitude", lwd=2,
          xlab = "Longitude", col = "black")
plot_ways(x = roads, add = TRUE, col = "grey")

points(bg_centroids, col = "magenta", type = "p", pch = 16, cex = 1)

box(lwd = 2)
grid()

## close output to file
dev.off()

#######################################
# write output
write.table(x = bg_centroids,
            file = paste0("./", gsub(pattern = " ", replacement = "_", x = title, fixed = TRUE), "_centroids.csv"),
            quote = FALSE, sep = ',', row.names = FALSE,
            col.names = c("Longitude", "Latitude"))












