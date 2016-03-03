### packages
require(dplyr)
require(ggplot2)
require(reshape2)
require(plyr)
require(lubridate)
require(sp)
require(rgdal)
require(maptools)
require(gridExtra)
require(RSEIS)
require(geosphere)
require(VTrack)
require(scales)

### creating pelagic dataset
wisefish_pel <- data.frame(depth=rep(NA, 70), size=rep(NA,70), dist=rep(NA, 70), distkm=rep(NA, 70),  distdeg=rep(NA, 70), shorelat=rep(NA,70), lat=rep(NA,70), long=rep(NA, 70))
wisefish_pel$dist[1:30] <- runif(n=30, min=3000, max=10000)
wisefish_pel$dist[31:50] <- runif(n=20, min=3000, max=5000)
wisefish_pel$dist[51:70] <- runif(n=20, min=3000, max=10000)
wisefish_pel$size <- runif(n=70, min=10, max=100)
wisefish_pel$depth <- wisefish_pel$dist * 0.02
wisefish_pel$distkm <- wisefish_pel$dist / 1000

wisefish_pel$long[1:30] <- runif(n=30, min=-64.29, max=-63.95)
wisefish_pel$long[31:50] <- runif(n=20, min=-63.95, max=-63.5)
wisefish_pel$long[51:70] <- runif(n=20, min=-64.3, max=-64.2)

wisefish_pel$shorelat[1:15] <- 45.385
wisefish_pel$shorelat[16:30] <- 45.25
wisefish_pel$shorelat[31:40] <- 45.385
wisefish_pel$shorelat[41:50] <- 45.315
wisefish_pel$shorelat[51:60] <- 45.15
wisefish_pel$shorelat[61:70] <- 45.275

wisefish_pel$distdeg <- wisefish_pel$distkm * 0.006

wisefish_pel$lat[1:15] <- wisefish_pel$shorelat[1:15] - wisefish_pel$distdeg[1:15]
wisefish_pel$lat[16:30] <- wisefish_pel$shorelat[16:30] + wisefish_pel$distdeg[16:30]
wisefish_pel$lat[31:40] <- wisefish_pel$shorelat[31:40] - wisefish_pel$distdeg[31:40]
wisefish_pel$lat[41:50] <- wisefish_pel$shorelat[41:50] + wisefish_pel$distdeg[41:50]
wisefish_pel$lat[51:60] <- wisefish_pel$shorelat[51:60] + wisefish_pel$distdeg[51:60]
wisefish_pel$lat[61:70] <- wisefish_pel$shorelat[61:70] - wisefish_pel$distdeg[61:70]

write.csv(wisefish_pel, "./wisefish_pel.csv")


### creating near shore dataset
wisefish_near <- data.frame(depth=rep(NA, 80), size=rep(NA,80), dist=rep(NA, 80), distkm=rep(NA, 80),  distdeg=rep(NA, 80), shorelat=rep(NA,80), shorelong=rep(NA, 80), lat=rep(NA,80), long=rep(NA, 80))
wisefish_near$dist <- runif(n=80, min=500, max=3000)
wisefish_near$size <- runif(n=80, min=5, max=60)
wisefish_near$depth <- wisefish_near$dist * 0.02
wisefish_near$distkm <- wisefish_near$dist / 1000
wisefish_near$distdeg <- wisefish_near$distkm * 0.006

shore_s <- subset(map, lat>45.12 & lat<45.35 & long > -64.39 & long < -63.5 & hole=="FALSE", select=c(long,lat))
shore_n <- subset(map, lat> 45.12 & lat>45.35 & long> -64.39 & long < -63.5 & hole=="FALSE", select=c(long,lat))
shore_s$unq <- c(1:348)
shore_n$unq <- c(1:1453)
fish_s <- data.frame(unq=rep(NA, 40))
fish_n <- data.frame(unq=rep(NA, 40))
fish_s$unq <- runif(n=40, min=1, max=348)
fish_n$unq <- runif(n=40, min=1, max=1453)
fish_s$unq <- round(fish_s$unq)
fish_n$unq <- round(fish_n$unq)
shore_s <- join(shore_s, fish_s, type="right")
shore_n <- join(shore_n, fish_n, type="right")

wisefish_near$shorelong[1:40] <- shore_s$long
wisefish_near$shorelat[1:40] <- shore_s$lat
wisefish_near$shorelong[41:80] <- shore_n$long
wisefish_near$shorelat[41:80] <- shore_n$lat

wisefish_near$lat[1:40] <- wisefish_near$distdeg[1:40] + wisefish_near$shorelat[1:40]
wisefish_near$lat[41:80] <- wisefish_near$shorelat[41:80] - wisefish_near$distdeg[41:80]
wisefish_near$long <- wisefish_near$shorelong

wisefish_near$long <- ifelse(wisefish_near$long< -64.3, wisefish_near$long +  wisefish_near$distdeg, wisefish_near$long)

write.csv(wisefish_near, "./wisefish_near.csv")

### reading in "good" datasets
wisefish_pel <- read.csv("./wisefish_pel.csv")
wisefish_near <- read.csv("./wisefish_near.csv")

wisefish_near <- select(wisefish_near, -shorelong)

wisefish_pel$type <- "pelagic"
wisefish_near$type <- "near shore"

wisefish <- rbind(wisefish_pel, wisefish_near)

### plotting on Minas Basin map
map_sp <- readShapeSpatial("./New_Polygon.shp")
map <- fortify(map_sp)
SP <- SpatialPoints(cbind(map$long, map$lat),
                    proj4string=CRS("+proj=utm +zone=20 +datum=WGS84 +units=m +no_defs "))
SP2 <- spTransform(SP, CRS("+proj=longlat"))

head(SP)
map <- cbind(map, SP2)
map$long <- map$coords.x1
map$lat <- map$coords.x2

xlim <- c(-64.8, -63.2)
ylim <- c(44.95, 45.52)


createScaleBar <- function(lon,lat,distanceLon,distanceLat,distanceLegend, dist.units = "km"){
  # First rectangle
  bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon, dist.units = dist.units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLat, dist.units = dist.units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"], bottomRight[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon*2, dist.units = dist.units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"], bottomRight2[1,"long"], bottomRight2[1,"long"], bottomRight[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLegend, dist.units = dist.units, model = "WGS84")
  onTop2 <- onTop3 <- onTop
  onTop2[1,"long"] <- bottomRight[1,"long"]
  onTop3[1,"long"] <- bottomRight2[1,"long"]
  
  legend <- rbind(onTop, onTop2, onTop3)
  legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}

createOrientationArrow <- function(scaleBar, length, distance = 1, dist.units = "km"){
  lon <- scaleBar$rectangle2[1,1]
  lat <- scaleBar$rectangle2[1,2]
  
  # Bottom point of the arrow
  begPoint <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist.units, model = "WGS84")
  lon <- begPoint[1,"long"]
  lat <- begPoint[1,"lat"]
  
  # Let us create the endpoint
  onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist.units, model = "WGS84")
  
  leftArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 225, dist = length/5, dist.units = dist.units, model = "WGS84")
  
  rightArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 135, dist = length/5, dist.units = dist.units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = onTop[1,"long"], yend = onTop[1,"lat"]),
    cbind(x = leftArrow[1,"long"], y = leftArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]),
    cbind(x = rightArrow[1,"long"], y = rightArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coordsN <- cbind(x = lon, y = (lat + onTop[1,"lat"])/2)
  
  return(list(res = res, coordsN = coordsN))
}

scaleBar <- function(lon, lat, distanceLon, distanceLat, distanceLegend, dist.unit = "km", rec.fill = "white", rec.colour = "black", rec2.fill = "black", rec2.colour = "black", legend.colour = "black", legend.size = 3, orientation = TRUE, arrow.length = 500, arrow.distance = 300, arrow.North.size = 6){
  laScaleBar <- createScaleBar(lon = lon, lat = lat, distanceLon = distanceLon, distanceLat = distanceLat, distanceLegend = distanceLegend, dist.unit = dist.unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = laScaleBar$rectangle, aes(x = lon, y = lat), fill = rec.fill, colour = rec.colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, aes(x = lon, y = lat), fill = rec2.fill, colour = rec2.colour)
  
  # Legend
  scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"], dist.unit, sep=""), x = laScaleBar$legend[,"long"], y = laScaleBar$legend[,"lat"], size = legend.size, colour = legend.colour)
  
  res <- list(rectangle1, rectangle2, scaleBarLegend)
  
  if(orientation){# Add an arrow pointing North
    coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, length = arrow.length, distance = arrow.distance, dist.unit = dist.unit)
    arrow <- list(geom_segment(data = coordsArrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coordsArrow$coordsN[1,"x"], y = coordsArrow$coordsN[1,"y"], size = arrow.North.size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}


ggplot() + geom_polygon(data=map, aes(x=long, y=lat, group=group, fill=hole)) + 
  geom_path(data=map, aes(x=long, y=lat, group=group, fill=hole), colour="black") +
  theme_bw(base_size=15) +
  coord_equal() +
  scale_fill_manual(values=c("lightgrey", "white"), guide="none") +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  theme(legend.background=element_rect(fill="white", colour="black", size=0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(colour="black")) +
  #geom_point(data=wisefish, aes(x=long, y=lat, shape=type))+
  xlab("Longitude") +
  ylab("Latitude") +
  #scale_shape_manual(values=c(4,1), guide=FALSE) +
  scale_x_continuous(breaks=c(-64.8, -64.6, -64.4, -64.2, -64.0, -63.8, -63.6, 
                              -63.4, -63.2))

wisefish <- select(wisefish, -X, -dist, -distdeg, -shorelat)

write.csv(wisefish, "./wisefish.csv")
wisefish <- read.csv("./wisefish.csv")
