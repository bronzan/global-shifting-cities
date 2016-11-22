 pythag <- function(lon1, lat1, lon2, lat2) {
 	x <- sqrt((lon2 - lon1)^2 + (lat2 - lat1)^2)
 	return(x)
 }


coord <- cities[1,2:3]

rowind <- which(pythag(coord$LON[1], coord$LAT[1], lonlat$lon, lonlat$lat) == min(pythag(coord$LON[1], coord$LAT[1], lonlat$lon, lonlat$lat)))

print(lonlat[rowind,])



 for (i in 1:nrow(marketcoord)) {
	rowind <- which(pythag(marketcoord[i,2], marketcoord[i,3], model_means$LON, model_means$LAT) == min(pythag(marketcoord[i,2], marketcoord[i,3], model_means$LON, model_means$LAT)))
	marketcoord$gridlon[i] <- model_means[rowind, 1]
	marketcoord$gridlat[i] <- model_means[rowind, 2]
}


newmarkets <- marketcoord
profiledata <- data.frame()


for (i in 1:nrow(newmarkets)) {
	stnlon <- marketcoord$gridlon[i]
	stnlat <- marketcoord$gridlat[i]
	newrow <- model_means[stnlon == model_means$LON & stnlat == model_means$LAT,3:367]
	profiledata <- rbind(profiledata, newrow)
}