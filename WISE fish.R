wisefish_pel <- data.frame(depth=rep(NA, 50), size=rep(NA,50), dist=rep(NA, 50), lat=rep(NA,50), long=rep(NA, 50))
wisefish_pel$depth <- rnorm(n=50, mean=50, sd=10)
wisefish_pel$size <- runif(n=50, min=10, max=100)
wisefish_pel$dist <- wisefish_pel$depth / 0.02
wisefish_pel$distkm <- wisefish_pel$dist / 1000

wisefish_pel$lat[1:25] <- runif(n=25, min=45.2, max=45.32)
wisefish_pel$lat[26:50] <- runif(n=25, min=45.32, max=45.4)
wisefish_pel$distlat <- wisefish_pel$distkm * 0.02
wisefish_pel$latloc[1:25] <- wisefish_pel$distlat[1:25] + wisefish_pel$lat[1:25]
wisefish_pel$latloc[26:50] <- wisefish_pel$lat[26:50] - wisefish_pel$distlat[26:50]

wisefish_pel <- arrange(wisefish_pel, dist)

wisefish_pel$long[1:5] <- runif(n=5, min=-63.65, max=-63.5)
wisefish_pel$long[6:10] <- runif(n=5, min=-63.85, max=-63.65)
wisefish_pel$long[11:50] <- runif(n=40, min=-64.3, max=-63.85)


wisefish_near <- data.frame(depth=rep(NA, 50), size=rep(NA,50), dist=rep(NA, 50), lat=rep(NA,50), long=rep(NA, 50))
wisefish_near$depth <- rnorm(n=50, mean=15, sd=5)
wisefish_near$size <- runif(n=50, min=5, max=60)
wisefish_near$dist <- wisefish_pel$depth / 0.02
wisefish_near$distkm <- wisefish_pel$dist / 1000

wisefish_near$lat[1:25] <- runif(n=25, min=45.2, max=45.32)
wisefish_near$lat[26:50] <- runif(n=25, min=45.32, max=45.4)
wisefish_near$distlat <- wisefish_near$distkm * 0.02
wisefish_near$latloc[1:25] <- wisefish_near$distlat[1:25] + wisefish_near$lat[1:25]
wisefish_near$latloc[26:50] <- wisefish_near$lat[26:50] - wisefish_near$distlat[26:50]

wisefish_near <- arrange(wisefish_near, dist)

wisefish_near$long[1:5] <- runif(n=5, min=-63.65, max=-63.5)
wisefish_near$long[6:10] <- runif(n=5, min=-63.85, max=-63.65)
wisefish_near$long[11:50] <- runif(n=40, min=-64.3, max=-63.85)




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



map_sp <- readShapeSpatial("/Volumes/Macintosh HD/Users/freyakeyser/Desktop/Shapefiles 2/New_Polygon.shp")
map <- fortify(map_sp)
SP <- SpatialPoints(cbind(map$long, map$lat),
                    proj4string=CRS("+proj=utm +zone=20 +datum=WGS84 +units=m +no_defs "))
SP2 <- spTransform(SP, CRS("+proj=longlat"))

head(SP)
map <- cbind(map, SP2)
map$long <- map$coords.x1
map$lat <- map$coords.x2

xlim <- c(-64.95, -63.2)
ylim <- c(44.95, 45.52)

ggplot() + geom_polygon(data=map, aes(x=long, y=lat, group=group, fill=hole)) + geom_path(data=map, aes(x=long, y=lat, group=group, fill=hole), colour="black") +
  theme_bw(base_size=15) +
  coord_equal() +
  scale_fill_manual(values=c("lightgrey", "white"), guide="none") +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  theme(legend.background=element_rect(fill="white", colour="black", size=0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(colour="black")) +
  geom_point(data=wisefish_pel, aes(x=long, y=latloc), shape=1, colour="blue")+
  xlab("Longitude") +
  ylab("Latitude")
