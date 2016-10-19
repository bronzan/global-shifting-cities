library(maps)
library(raster)
library(geosphere)

# FROM CMIP5

# for Southern Hemisphere summers
cities.now.2.6.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg26DJFhist.nc"
cities.then.2.6.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg26DJF.nc"

cities.now.4.5.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg45DJFhist.nc"
cities.then.4.5.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg45DJF.nc"

cities.now.6.0.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg60DJFhist.nc"
cities.then.6.0.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg60DJF.nc"

cities.now.8.5.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg85DJFhist.nc"
cities.then.8.5.DJF.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg85DJF.nc"

r.now.2.6.DJF <- raster(cities.now.2.6.DJF.file)
r.then.2.6.DJF <- raster(cities.then.2.6.DJF.file)
r.diff.2.6.DJF <- r.then.2.6.DJF - r.now.2.6.DJF

r.now.4.5.DJF <- raster(cities.now.4.5.DJF.file)
r.then.4.5.DJF <- raster(cities.then.4.5.DJF.file)
r.diff.4.5.DJF <- r.then.4.5.DJF - r.now.4.5.DJF

r.now.6.0.DJF <- raster(cities.now.6.0.DJF.file)
r.then.6.0.DJF <- raster(cities.then.6.0.DJF.file)
r.diff.6.0.DJF <- r.then.6.0.DJF - r.now.6.0.DJF

r.now.8.5.DJF <- raster(cities.now.8.5.DJF.file)
r.then.8.5.DJF <- raster(cities.then.8.5.DJF.file)
r.diff.8.5.DJF <- r.then.8.5.DJF - r.now.8.5.DJF

# for Northern Hemisphere summers

cities.now.2.6.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg26JJAhist.nc"
cities.then.2.6.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg26JJA.nc"

cities.now.4.5.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg45JJAhist.nc"
cities.then.4.5.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg45JJA.nc"

cities.now.6.0.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg60JJAhist.nc"
cities.then.6.0.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg60JJA.nc"

cities.now.8.5.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg85JJAhist.nc"
cities.then.8.5.JJA.file <- "/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/avg85JJA.nc"

r.now.2.6.JJA <- raster(cities.now.2.6.JJA.file)
r.then.2.6.JJA <- raster(cities.then.2.6.JJA.file)
r.diff.2.6.JJA <- r.then.2.6.JJA - r.now.2.6.JJA

r.now.4.5.JJA <- raster(cities.now.4.5.JJA.file)
r.then.4.5.JJA <- raster(cities.then.4.5.JJA.file)
r.diff.4.5.JJA <- r.then.4.5.JJA - r.now.4.5.JJA

r.now.6.0.JJA <- raster(cities.now.6.0.JJA.file)
r.then.6.0.JJA <- raster(cities.then.6.0.JJA.file)
r.diff.6.0.JJA <- r.then.6.0.JJA - r.now.6.0.JJA

r.now.8.5.JJA <- raster(cities.now.8.5.JJA.file)
r.then.8.5.JJA <- raster(cities.then.8.5.JJA.file)
r.diff.8.5.JJA <- r.then.8.5.JJA - r.now.8.5.JJA


shift.cities <- read.csv("~/Desktop/Tables/globalSummerMeans.csv")
shift.cities$NAME <- trimws(shift.cities$NAME)
shift.cities$elong <- shift.cities$LONGITUDE
shift.cities$elong[shift.cities$LONGITUDE<0]<-360+shift.cities$LONGITUDE[shift.cities$LONGITUDE<0]
country.codes <- read.csv("~/Desktop/Tables/country-codes.csv")
shift.cities <- merge(shift.cities, country.codes, by.x="countryCode", by.y="code", all.x=TRUE)

# 2a.  find model anomalies in CMIP5 JJA data for coordinates in northern hemisphere
shift.cities$diff.2.6[shift.cities$summerSzn=='JJA'] <- extract(r.diff.2.6.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.4.5[shift.cities$summerSzn=='JJA'] <- extract(r.diff.4.5.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.6.0[shift.cities$summerSzn=='JJA'] <- extract(r.diff.6.0.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.8.5[shift.cities$summerSzn=='JJA'] <- extract(r.diff.8.5.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')

# 2b.  find model anomalies in CMIP5 DJF data for coordinates in southern hemisphere
shift.cities$diff.2.6[shift.cities$summerSzn=='DJF'] <- extract(r.diff.2.6.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.4.5[shift.cities$summerSzn=='DJF'] <- extract(r.diff.4.5.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.6.0[shift.cities$summerSzn=='DJF'] <- extract(r.diff.6.0.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

shift.cities$diff.8.5[shift.cities$summerSzn=='DJF'] <- extract(r.diff.8.5.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')


# ALTERNATE EXTRACT FOR COASTAL NA'S

shift.cities$diff.2.6.nafix[shift.cities$summerSzn=='JJA'] <- extract(r.diff.2.6.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA', c("elong","LATITUDE")]), fun=mean, buffer=50000)
shift.cities$diff.2.6.nafix[shift.cities$summerSzn=='DJF'] <- extract(r.diff.2.6.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF', c("elong","LATITUDE")]), fun=mean, buffer=50000)

shift.cities$diff.4.5.nafix[shift.cities$summerSzn=='JJA'] <- extract(r.diff.4.5.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA', c("elong","LATITUDE")]), fun=mean, buffer=50000)
shift.cities$diff.4.5.nafix[shift.cities$summerSzn=='DJF'] <- extract(r.diff.4.5.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF', c("elong","LATITUDE")]), fun=mean, buffer=50000)

shift.cities$diff.6.0.nafix[shift.cities$summerSzn=='JJA'] <- extract(r.diff.6.0.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA', c("elong","LATITUDE")]), fun=mean, buffer=50000)
shift.cities$diff.6.0.nafix[shift.cities$summerSzn=='DJF'] <- extract(r.diff.6.0.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF', c("elong","LATITUDE")]), fun=mean, buffer=50000)

shift.cities$diff.8.5.nafix[shift.cities$summerSzn=='JJA'] <- extract(r.diff.8.5.JJA,  SpatialPoints(shift.cities[shift.cities$summerSzn=='JJA', c("elong","LATITUDE")]), fun=mean, buffer=50000)
shift.cities$diff.8.5.nafix[shift.cities$summerSzn=='DJF'] <- extract(r.diff.8.5.DJF,  SpatialPoints(shift.cities[shift.cities$summerSzn=='DJF', c("elong","LATITUDE")]), fun=mean, buffer=50000)

shift.cities$diff.2.6[is.na(shift.cities$diff.2.6)] <- shift.cities$diff.2.6.nafix[is.na(shift.cities$diff.2.6)]
shift.cities$diff.4.5[is.na(shift.cities$diff.4.5)] <- shift.cities$diff.4.5.nafix[is.na(shift.cities$diff.4.5)]
shift.cities$diff.6.0[is.na(shift.cities$diff.6.0)] <- shift.cities$diff.6.0.nafix[is.na(shift.cities$diff.6.0)]
shift.cities$diff.8.5[is.na(shift.cities$diff.8.5)] <- shift.cities$diff.8.5.nafix[is.na(shift.cities$diff.8.5)]


shift.cities$now <- shift.cities$summerMean

# 3.  add anomolies to temps
shift.cities$then.2.6 <- shift.cities$now + shift.cities$diff.2.6

shift.cities$then.4.5 <- shift.cities$now + shift.cities$diff.4.5

shift.cities$then.6.0 <- shift.cities$now + shift.cities$diff.6.0

shift.cities$then.8.5 <- shift.cities$now + shift.cities$diff.8.5

# Luxor, Egypt: 41.1 (on par w/ Lake Havasu, AZ)
# Khartoum, Sudan: 41.7
# Abu Dhabi, UAE: 43
# Riyadh, Saudi Arabia: 43.9
# And...Kuwait City, Kuwait: 45.6

mideast <- data.frame(rbind(c("Luxor", "Egypt", 25.70, 32.65, 41.1),
				c("Khartoum", "Sudan", 15.58, 32.52, 41.7),
				c("Abu Dhabi", "UAE", 24.48, 54.37, 43),
				c("Riyadh", "Saudi Arabia", 24.65, 46.77, 43.9),
				c("Kuwait City", "Kuwait", 29.38, 47.99, 45.6)))
colnames(mideast) <- c("name", "country", "lat", "long", "now")
class(mideast$now) <- "numeric"

me.cities <- c("Luxor", "Khartoum", "Abu Dhabi", "Riyadh", "Kuwait City")
# subset(world.cities, name %in% me.cities)
# 				           name          country.etc     pop   lat  long capital
# 265   Abu Dhabi United Arab Emirates  619316 24.48 54.37       1
# 17883  Khartoum                Sudan 2090001 15.58 32.52       1
# 21793     Luxor                Egypt  430214 25.70 32.65       0
# 31409    Riyadh         Saudi Arabia 4328067 24.65 46.77       1
# subset(world.cities, country.etc == "Kuwait" & capital == 1)
#           name country.etc   pop   lat  long capital
# 43361 al-Kuwayt      Kuwait 63596 29.38 47.99       1


match.cities <- shift.cities

all.city.matches.out <- data.frame()
match.3step <- data.frame()
match.hemi <- data.frame()
straight.temp.match <- data.frame()
# 4.  capitals loop
# for (i in 1:nrow(us.capitals)){
for (i in 1:nrow(shift.cities)){

	# capital <- us.capitals[i,]
	capital <- shift.cities[i,]

# 	a.  find other cities that are SSE-S-SSW of current city
	cap.bearing <- bearing(capital[,c("LONGITUDE","LATITUDE")], match.cities[,c("LONGITUDE","LATITUDE")])

# 	b.  calculate distances
#	cap.distance <- distHaversine(capital[,c("long","lat")], match.cities[,c("long","lat")], r=6378137)

#	b.2. calculate latitude and longitude differences
	cap.lat.diff <- abs(abs(capital$LATITUDE) - abs(match.cities$LATITUDE))
	cap.lat.diff.h <- cap.lat.diff
	cap.lat.diff.h[(match.cities$LATITUDE < 0) != (capital$LATITUDE < 0)] <- 0

#	cap.long.diff <- capital$LONGITUDE - match.cities$LONGITUDE

# 	c.  calculate temperature difference
	cap.diff.2.6 <- capital$then.2.6 - match.cities$now
	cap.diff.4.5 <- capital$then.4.5 - match.cities$now
	cap.diff.6.0 <- capital$then.6.0 - match.cities$now
	cap.diff.8.5 <- capital$then.8.5 - match.cities$now

#	d.  select all with temp match within 0.5 degrees then look at "best" bearing

	temp.order.cities.2.6 <- match.cities[order(abs(cap.diff.2.6), decreasing=FALSE),]
	temp.order.cities.4.5 <- match.cities[order(abs(cap.diff.4.5), decreasing=FALSE),]
	temp.order.cities.6.0 <- match.cities[order(abs(cap.diff.6.0), decreasing=FALSE),]
	temp.order.cities.8.5 <- match.cities[order(abs(cap.diff.8.5), decreasing=FALSE),]

	temp.order.bearing.2.6 <- cap.bearing[order(abs(cap.diff.2.6), decreasing=FALSE)]
	temp.order.bearing.4.5 <- cap.bearing[order(abs(cap.diff.4.5), decreasing=FALSE)]
	temp.order.bearing.6.0 <- cap.bearing[order(abs(cap.diff.6.0), decreasing=FALSE)]
	temp.order.bearing.8.5 <- cap.bearing[order(abs(cap.diff.8.5), decreasing=FALSE)]

	temp.order.lat.diff.2.6 <- cap.lat.diff[order(abs(cap.diff.2.6), decreasing=FALSE)]
	temp.order.lat.diff.4.5 <- cap.lat.diff[order(abs(cap.diff.4.5), decreasing=FALSE)]
	temp.order.lat.diff.6.0 <- cap.lat.diff[order(abs(cap.diff.6.0), decreasing=FALSE)]
	temp.order.lat.diff.8.5 <- cap.lat.diff[order(abs(cap.diff.8.5), decreasing=FALSE)]

	temp.order.lat.diff.h.2.6 <- cap.lat.diff.h[order(abs(cap.diff.2.6), decreasing=FALSE)]
	temp.order.lat.diff.h.4.5 <- cap.lat.diff.h[order(abs(cap.diff.4.5), decreasing=FALSE)]
	temp.order.lat.diff.h.6.0 <- cap.lat.diff.h[order(abs(cap.diff.6.0), decreasing=FALSE)]
	temp.order.lat.diff.h.8.5 <- cap.lat.diff.h[order(abs(cap.diff.8.5), decreasing=FALSE)]

	if (capital$LATITUDE <= 0){
		# small bearings will be "north" for SH cities
		closest.temp.bearing.2.6 <- abs(temp.order.bearing.2.6[1:10])
		closest.temp.bearing.4.5 <- abs(temp.order.bearing.4.5[1:10])
		closest.temp.bearing.6.0 <- abs(temp.order.bearing.6.0[1:10])
		closest.temp.bearing.8.5 <- abs(temp.order.bearing.8.5[1:10])
	} else {
		# small bearings will be "south" for NH cities
		closest.temp.bearing.2.6 <- abs( abs(temp.order.bearing.2.6[1:10]) - 180 )
		closest.temp.bearing.4.5 <- abs( abs(temp.order.bearing.4.5[1:10]) - 180 )
		closest.temp.bearing.6.0 <- abs( abs(temp.order.bearing.6.0[1:10]) - 180 )
		closest.temp.bearing.8.5 <- abs( abs(temp.order.bearing.8.5[1:10]) - 180 )
	}
### NEW : TOP 10 CLOSEST BY TEMP
###		  TOP 3 BEST BEARING (toward equator)
###		  TOP 1 MOST SOUTH (biggest lat diff)

	bearing.three.2.6 <- temp.order.cities.2.6[1:10,][order(closest.temp.bearing.2.6, decreasing=FALSE),][1:3,]
	bearing.three.4.5 <- temp.order.cities.4.5[1:10,][order(closest.temp.bearing.4.5, decreasing=FALSE),][1:3,]
	bearing.three.6.0 <- temp.order.cities.6.0[1:10,][order(closest.temp.bearing.6.0, decreasing=FALSE),][1:3,]
	bearing.three.8.5 <- temp.order.cities.8.5[1:10,][order(closest.temp.bearing.8.5, decreasing=FALSE),][1:3,]

	lat.diff.three.2.6 <- temp.order.lat.diff.2.6[1:10][order(closest.temp.bearing.2.6, decreasing=FALSE)][1:3]
	lat.diff.three.4.5 <- temp.order.lat.diff.4.5[1:10][order(closest.temp.bearing.4.5, decreasing=FALSE)][1:3]
	lat.diff.three.6.0 <- temp.order.lat.diff.6.0[1:10][order(closest.temp.bearing.6.0, decreasing=FALSE)][1:3]
	lat.diff.three.8.5 <- temp.order.lat.diff.8.5[1:10][order(closest.temp.bearing.8.5, decreasing=FALSE)][1:3]

	lat.diff.h.three.2.6 <- temp.order.lat.diff.h.2.6[1:10][order(closest.temp.bearing.2.6, decreasing=FALSE)][1:3]
	lat.diff.h.three.4.5 <- temp.order.lat.diff.h.4.5[1:10][order(closest.temp.bearing.4.5, decreasing=FALSE)][1:3]
	lat.diff.h.three.6.0 <- temp.order.lat.diff.h.6.0[1:10][order(closest.temp.bearing.6.0, decreasing=FALSE)][1:3]
	lat.diff.h.three.8.5 <- temp.order.lat.diff.h.8.5[1:10][order(closest.temp.bearing.8.5, decreasing=FALSE)][1:3]

	mtch.2.6 <- bearing.three.2.6[which.max(lat.diff.three.2.6),]
	mtch.4.5 <- bearing.three.4.5[which.max(lat.diff.three.4.5),]
	mtch.6.0 <- bearing.three.6.0[which.max(lat.diff.three.6.0),]
	mtch.8.5 <- bearing.three.8.5[which.max(lat.diff.three.8.5),]

	mtch.h.2.6 <- bearing.three.2.6[which.max(lat.diff.h.three.2.6),]
	mtch.h.4.5 <- bearing.three.4.5[which.max(lat.diff.h.three.4.5),]
	mtch.h.6.0 <- bearing.three.6.0[which.max(lat.diff.h.three.6.0),]
	mtch.h.8.5 <- bearing.three.8.5[which.max(lat.diff.h.three.8.5),]

	straight.2.6 <- temp.order.cities.2.6[1,]
	straight.4.5 <- temp.order.cities.4.5[1,]
	straight.6.0 <- temp.order.cities.6.0[1,]
	straight.8.5 <- temp.order.cities.8.5[1,]

	print(paste0(capital$NAME, ", ", capital$LATITUDE, " - ", mtch.2.6$NAME, ", ", mtch.2.6$LATITUDE, " - ", mtch.4.5$NAME, ", ", mtch.4.5$LATITUDE, " - ", mtch.6.0$NAME, ", ", mtch.6.0$LATITUDE, " - ", mtch.8.5$NAME, ", ", mtch.8.5$LATITUDE))

	# print(paste0(capital$name, ", ", capital$country.etc, " - ", dim(ordered.cities.2.6)[1], " - ", dim(ordered.cities.4.5)[1], " - ", dim(ordered.cities.6.0)[1], dim(ordered.cities.8.5)[1]))

	all.city.matches.out <- rbind(all.city.matches.out, c(city=capital$name, mtch.2.6=mtch.2.6$name, mtch.4.5=mtch.4.5$name, mtch.6.0=mtch.6.0$name, mtch.8.5=mtch.8.5$name))

	one.3step.row <- c(capital$NAME, capital$country, capital$LATITUDE, capital$LONGITUDE, round(capital$now,2), mtch.2.6$NAME, mtch.2.6$country, mtch.2.6$LATITUDE, mtch.2.6$LONGITUDE, round(mtch.2.6$now,2), round(capital$then.2.6,2), mtch.4.5$NAME, mtch.4.5$country, mtch.4.5$LATITUDE, mtch.4.5$LONGITUDE, round(mtch.4.5$now,2), round(capital$then.4.5,2), mtch.6.0$NAME, mtch.6.0$country, mtch.6.0$LATITUDE, mtch.6.0$LONGITUDE, round(mtch.6.0$now,2), round(capital$then.6.0,2), mtch.8.5$NAME, mtch.8.5$country, mtch.8.5$LATITUDE, mtch.8.5$LONGITUDE, round(mtch.8.5$now,2), round(capital$then.8.5,2))

	one.hemi.row <- c(capital$NAME, capital$country, capital$LATITUDE, capital$LONGITUDE, round(capital$now,2), mtch.h.2.6$NAME, mtch.h.2.6$country, mtch.h.2.6$LATITUDE, mtch.h.2.6$LONGITUDE, round(mtch.h.2.6$now,2), round(capital$then.2.6,2), mtch.h.4.5$NAME, mtch.h.4.5$country, mtch.h.4.5$LATITUDE, mtch.h.4.5$LONGITUDE, round(mtch.h.4.5$now,2), round(capital$then.4.5,2), mtch.h.6.0$NAME, mtch.h.6.0$country, mtch.h.6.0$LATITUDE, mtch.h.6.0$LONGITUDE, round(mtch.h.6.0$now,2), round(capital$then.6.0,2), mtch.h.8.5$NAME, mtch.h.8.5$country, mtch.h.8.5$LATITUDE, mtch.h.8.5$LONGITUDE, round(mtch.h.8.5$now,2), round(capital$then.8.5,2))

	one.straight.row <- c(capital$NAME, capital$country, capital$LATITUDE, capital$LONGITUDE, round(capital$now,2), straight.2.6$NAME, straight.2.6$country, straight.2.6$LATITUDE, straight.2.6$LONGITUDE, round(straight.2.6$now,2), round(capital$then.2.6,2), straight.4.5$NAME, straight.4.5$country, straight.4.5$LATITUDE, straight.4.5$LONGITUDE, round(straight.4.5$now,2), round(capital$then.4.5,2), straight.6.0$NAME, straight.6.0$country, straight.6.0$LATITUDE, straight.6.0$LONGITUDE, round(straight.6.0$now,2), round(capital$then.6.0,2), straight.8.5$NAME, straight.8.5$country, straight.8.5$LATITUDE, straight.8.5$LONGITUDE, round(straight.8.5$now,2), round(capital$then.8.5,2))

	if (!is.na(capital$then.8.5) & capital$then.8.5 > 41.48413) {
		me.diff <- abs(capital$then.8.5 - mideast$now)
		me.match <- mideast[which.min(me.diff),]

		one.3step.row <- c(capital$NAME, capital$country, capital$LATITUDE, capital$LONGITUDE, round(capital$now,2), mtch.2.6$NAME, mtch.2.6$country, mtch.2.6$LATITUDE, mtch.2.6$LONGITUDE, round(mtch.2.6$now,2), round(capital$then.2.6,2), mtch.4.5$NAME, mtch.4.5$country, mtch.4.5$LATITUDE, mtch.4.5$LONGITUDE, round(mtch.4.5$now,2), round(capital$then.4.5,2), mtch.6.0$NAME, mtch.6.0$country, mtch.6.0$LATITUDE, mtch.6.0$LONGITUDE, round(mtch.6.0$now,2), round(capital$then.6.0,2), me.match$name, me.match$country, me.match$lat, me.match$long, round(me.match$now,2), round(capital$then.8.5,2))
	}
	match.3step <- rbind(match.3step, one.3step.row)
	match.hemi <- rbind(match.hemi, one.hemi.row)
	straight.temp.match <- rbind(straight.temp.match, one.straight.row)
}

colnames(match.3step) <- c("Station.Name", "Country", "Latitude", "Longitude", "Temperature", "Compared.2.6.Station.Name", "Compared.2.6.Country", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.Station.Name", "Compared.4.5.Country", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.Station.Name", "Compared.6.0.Country", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.Station.Name", "Compared.8.5.Country", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature")
colnames(match.hemi) <- c("Station.Name", "Country", "Latitude", "Longitude", "Temperature", "Compared.2.6.Station.Name", "Compared.2.6.Country", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.Station.Name", "Compared.4.5.Country", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.Station.Name", "Compared.6.0.Country", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.Station.Name", "Compared.8.5.Country", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature")
colnames(straight.temp.match) <- c("Station.Name", "Country", "Latitude", "Longitude", "Temperature", "Compared.2.6.Station.Name", "Compared.2.6.Country", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.Station.Name", "Compared.4.5.Country", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.Station.Name", "Compared.6.0.Country", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.Station.Name", "Compared.8.5.Country", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature")

write.csv(match.3step, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_3stepMatch.csv", row.names=FALSE)
write.csv(match.hemi, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_HemiMatch.csv", row.names=FALSE)
write.csv(straight.temp.match, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_StraightTempMatch.csv", row.names=FALSE)


lower.2.6 <- min(shift.cities$then.2.6 - 0.5, na.rm=TRUE)
upper.2.6 <- max(shift.cities$then.2.6 + 0.5, na.rm=TRUE)
buckets.2.6 <- as.numeric(cut(shift.cities$now[shift.cities$now > lower.2.6 & shift.cities$now < upper.2.6], 20))

lower.4.5 <- min(shift.cities$then.4.5 - 0.5, na.rm=TRUE)
upper.4.5 <- max(shift.cities$then.4.5 + 0.5, na.rm=TRUE)
buckets.4.5 <- as.numeric(cut(shift.cities$now[shift.cities$now > lower.4.5 & shift.cities$now < upper.4.5], 20))

lower.6.0 <- min(shift.cities$then.6.0 - 0.5, na.rm=TRUE)
upper.6.0 <- max(shift.cities$then.6.0 + 0.5, na.rm=TRUE)
buckets.6.0 <- as.numeric(cut(shift.cities$now[shift.cities$now > lower.6.0 & shift.cities$now < upper.6.0], 20))

lower.8.5 <- min(shift.cities$then.8.5 - 0.5, na.rm=TRUE)
upper.8.5 <- max(shift.cities$then.8.5 + 0.5, na.rm=TRUE)
buckets.8.5 <- as.numeric(cut(shift.cities$now[shift.cities$now > lower.8.5 & shift.cities$now < upper.8.5], 20))

cbind(buckets.2.6, buckets.4.5, buckets.6.0, buckets.8.5)

buckets <- cut(shift.cities$now, 40, dig.lab = 4)
bucket.levels <- as.numeric(buckets)
bucket.cities <- cbind(shift.cities, buckets, bucket.levels)

bucket.labs <- levels(buckets)
bucket.points <- data.frame(lower = as.numeric( sub("\\((.+),.*", "\\1", bucket.labs) ), upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", bucket.labs) ))

##########################################################
##  HALF-DEGREE BUCKETS - ROUND TO NEAREST 1/2 DEGREE?  ##
##########################################################

# bucket.cities$halfDegree <- round(bucket.cities$now*2)/2

rat.cities <- read.csv("~/Desktop/Tables/The Upshot - global shifting cities - buckets.rationalized.csv")

rat2.cities <- merge(rat.cities[,c("NAME", "country", "Actual.City", "countryCode", "WMOID", "ID", "destination")], shift.cities, all.x=TRUE)

rat.cities <- rat2.cities

rat.cities$half.now <- round(rat.cities$now*2)/2
rat.cities$half.2.6 <- round(rat.cities$then.2.6*2)/2
rat.cities$half.4.5 <- round(rat.cities$then.4.5*2)/2
rat.cities$half.6.0 <- round(rat.cities$then.6.0*2)/2
rat.cities$half.8.5 <- round(rat.cities$then.8.5*2)/2

rat.cities$floor.2.6 <- floor(rat.cities$then.2.6*2)/2
rat.cities$floor.4.5 <- floor(rat.cities$then.4.5*2)/2
rat.cities$floor.6.0 <- floor(rat.cities$then.6.0*2)/2
rat.cities$floor.8.5 <- floor(rat.cities$then.8.5*2)/2

rat.cities$ceil.2.6 <- ceil(rat.cities$then.2.6*2)/2
rat.cities$ceil.4.5 <- ceil(rat.cities$then.4.5*2)/2
rat.cities$ceil.6.0 <- ceil(rat.cities$then.6.0*2)/2
rat.cities$ceil.8.5 <- ceil(rat.cities$then.8.5*2)/2

rat.cities$alt.2.6 <- ifelse(rat.cities$half.2.6 == rat.cities$floor.2.6, rat.cities$ceil.2.6, rat.cities$floor.2.6 )
rat.cities$alt.4.5 <- ifelse(rat.cities$half.4.5 == rat.cities$floor.4.5, rat.cities$ceil.4.5, rat.cities$floor.4.5 )
rat.cities$alt.6.0 <- ifelse(rat.cities$half.6.0 == rat.cities$floor.6.0, rat.cities$ceil.6.0, rat.cities$floor.6.0 )
rat.cities$alt.8.5 <- ifelse(rat.cities$half.8.5 == rat.cities$floor.8.5, rat.cities$ceil.8.5, rat.cities$floor.8.5 )


rat.dest <- subset(rat.cities, bucket.destination==1)
# rat.dest <- rat.cities[!duplicated(rat.cities$half.now),]

bucket.match <- data.frame()
for (i in 1:nrow(rat.cities)){
	rat.city <- rat.cities[i,]

	bucket.2.6 <- rat.dest[rat.city$half.2.6 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
	if (nrow(bucket.2.6) == 0) {
		bucket.2.6 <- rat.dest[rat.city$half.2.6 == rat.dest$half.now, ]
		if (nrow(bucket.2.6) == 0) {
			bucket.2.6 <- rat.dest[rat.city$alt.2.6 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
			if (nrow(bucket.2.6) == 0) {
				bucket.2.6 <- rat.dest[rat.city$alt.2.6 == rat.dest$half.now, ]
			}
		}
	}

	bucket.4.5 <- rat.dest[rat.city$half.4.5 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
	if (nrow(bucket.4.5) == 0) {
		bucket.4.5 <- rat.dest[rat.city$half.4.5 == rat.dest$half.now, ]
		if (nrow(bucket.4.5) == 0) {
			bucket.4.5 <- rat.dest[rat.city$alt.4.5 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
			if (nrow(bucket.4.5) == 0) {
				bucket.4.5 <- rat.dest[rat.city$alt.4.5 == rat.dest$half.now, ]
			}
		}
	}

	bucket.6.0 <- rat.dest[rat.city$half.6.0 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
	if (nrow(bucket.6.0) == 0) {
		bucket.6.0 <- rat.dest[rat.city$half.6.0 == rat.dest$half.now, ]
		if (nrow(bucket.6.0) == 0) {
			bucket.6.0 <- rat.dest[rat.city$alt.6.0 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
			if (nrow(bucket.6.0) == 0) {
				bucket.6.0 <- rat.dest[rat.city$alt.6.0 == rat.dest$half.now, ]
			}
		}
	}

	bucket.8.5 <- rat.dest[rat.city$half.8.5 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
	if (nrow(bucket.8.5) == 0) {
		bucket.8.5 <- rat.dest[rat.city$half.8.5 == rat.dest$half.now, ]
		if (nrow(bucket.8.5) == 0) {
			bucket.8.5 <- rat.dest[rat.city$alt.8.5 == rat.dest$half.now & rat.city$summerSzn == rat.dest$summerSzn, ]
			if (nrow(bucket.8.5) == 0) {
				bucket.8.5 <- rat.dest[rat.city$alt.8.5 == rat.dest$half.now, ]
			}
		}
	}


	one.bucket.row <- c(rat.city$Actual.City, rat.city$country, rat.city$LATITUDE,  rat.city$LONGITUDE, round(rat.city$now,2), ifelse(length(bucket.2.6$Actual.City) > 0, bucket.2.6$Actual.City, NA), ifelse(length(bucket.2.6$country) > 0, bucket.2.6$country, NA), ifelse(length(bucket.2.6$LATITUDE) > 0, bucket.2.6$LATITUDE, NA), ifelse(length(bucket.2.6$LONGITUDE) > 0, bucket.2.6$LONGITUDE, NA), ifelse(length(round(bucket.2.6$now,2)) > 0, round(bucket.2.6$now,2), NA), ifelse(length(round(rat.city$then.2.6,2)) > 0, round(rat.city$then.2.6,2), NA), ifelse(length(bucket.4.5$Actual.City) > 0, bucket.4.5$Actual.City, NA), ifelse(length(bucket.4.5$country) > 0, bucket.4.5$country, NA), ifelse(length(bucket.4.5$LATITUDE) > 0, bucket.4.5$LATITUDE, NA), ifelse(length(bucket.4.5$LONGITUDE) > 0, bucket.4.5$LONGITUDE, NA), ifelse(length(round(bucket.4.5$now,2)) > 0, round(bucket.4.5$now,2), NA), ifelse(length(round(rat.city$then.4.5,2)) > 0, round(rat.city$then.4.5,2), NA), ifelse(length(bucket.6.0$Actual.City) > 0, bucket.6.0$Actual.City, NA), ifelse(length(bucket.6.0$country) > 0, bucket.6.0$country, NA), ifelse(length(bucket.6.0$LATITUDE) > 0, bucket.6.0$LATITUDE, NA), ifelse(length(bucket.6.0$LONGITUDE) > 0, bucket.6.0$LONGITUDE, NA), ifelse(length(round(bucket.6.0$now,2)) > 0, round(bucket.6.0$now,2), NA), ifelse(length(round(rat.city$then.6.0,2)) > 0, round(rat.city$then.6.0,2), NA), ifelse(length(bucket.8.5$Actual.City) > 0, bucket.8.5$Actual.City, NA), ifelse(length(bucket.8.5$country) > 0, bucket.8.5$country, NA), ifelse(length(bucket.8.5$LATITUDE) > 0, bucket.8.5$LATITUDE, NA), ifelse(length(bucket.8.5$LONGITUDE) > 0, bucket.8.5$LONGITUDE, NA), ifelse(length(round(bucket.8.5$now,2)) > 0, round(bucket.8.5$now,2), NA), ifelse(length(round(rat.city$then.8.5,2)) > 0, round(rat.city$then.8.5,2), NA) )

	bucket.match <- rbind(bucket.match, one.bucket.row)
}

colnames(bucket.match) <- c( "City.Name", "Country", "Latitude", "Longitude", "Temperature",
	"Compared.2.6.City.Name", "Compared.2.6.Country", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.City.Name", "Compared.4.5.Country", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature",
 	"Compared.6.0.City.Name", "Compared.6.0.Country", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.City.Name", "Compared.8.5.Country", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature" )

write.csv(bucket.match, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_bucketMatch.csv", row.names=FALSE)




### USE BUCKET "DESTINATIONS" TO DO AN EXACT MATCH


rat.dest <- subset(rat.cities, destination==1)

straight.dest.match <- data.frame()
for (i in 1:nrow(rat.cities)){
	rat.city <- rat.cities[i,]

	bucket.2.6 <- rat.dest[which.min(abs(rat.city$then.2.6 - rat.dest$now)), ]
	if (nrow(bucket.2.6) > 1)
	{
		bucket.2.6.temp <- bucket.2.6[bucket.2.6$summerSzn == rat.city$summzerSzn]
		if (nrow(bucket.2.6.temp) == 0) {bucket.2.6.temp <- bucket.2.6}
		bucket.2.6 <- bucket.2.6.temp
		if (nrow(bucket.2.6) > 1) {bucket.2.6.temp <- bucket.2.6[which.min(abs(bucket.2.6$LONGITUDE - rat.city$LONGITUDE)),]}
	}

	bucket.4.5 <- rat.dest[which.min(abs(rat.city$then.4.5 - rat.dest$now)), ]
	if (nrow(bucket.4.5) > 1)
	{
		bucket.4.5.temp <- bucket.4.5[bucket.4.5$summerSzn == rat.city$summzerSzn]
		if (nrow(bucket.4.5.temp) == 0) {bucket.4.5.temp <- bucket.4.5}
		bucket.4.5 <- bucket.4.5.temp
		if (nrow(bucket.4.5) > 1) {bucket.4.5.temp <- bucket.4.5[which.min(abs(bucket.4.5$LONGITUDE - rat.city$LONGITUDE)),]}
	}

	bucket.6.0 <- rat.dest[which.min(abs(rat.city$then.6.0 - rat.dest$now)), ]
	if (nrow(bucket.6.0) > 1)
	{
		bucket.6.0.temp <- bucket.6.0[bucket.6.0$summerSzn == rat.city$summzerSzn]
		if (nrow(bucket.6.0.temp) == 0) {bucket.6.0.temp <- bucket.6.0}
		bucket.6.0 <- bucket.6.0.temp
		if (nrow(bucket.6.0) > 1) {bucket.6.0.temp <- bucket.6.0[which.min(abs(bucket.6.0$LONGITUDE - rat.city$LONGITUDE)),]}
	}

	bucket.8.5 <- rat.dest[which.min(abs(rat.city$then.8.5 - rat.dest$now)), ]
	if (nrow(bucket.8.5) > 1)
	{
		bucket.8.5.temp <- bucket.8.5[bucket.8.5$summerSzn == rat.city$summzerSzn]
		if (nrow(bucket.8.5.temp) == 0) {bucket.8.5.temp <- bucket.8.5}
		bucket.8.5 <- bucket.8.5.temp
		if (nrow(bucket.8.5) > 1) {bucket.8.5.temp <- bucket.8.5[which.min(abs(bucket.8.5$LONGITUDE - rat.city$LONGITUDE)),]}
	}


	one.dest.row <- c(rat.city$Actual.City, rat.city$country, rat.city$LATITUDE,  rat.city$LONGITUDE, round(rat.city$now,2), ifelse(length(bucket.2.6$Actual.City) > 0, bucket.2.6$Actual.City, NA), ifelse(length(bucket.2.6$country) > 0, bucket.2.6$country, NA), ifelse(length(bucket.2.6$LATITUDE) > 0, bucket.2.6$LATITUDE, NA), ifelse(length(bucket.2.6$LONGITUDE) > 0, bucket.2.6$LONGITUDE, NA), ifelse(length(round(bucket.2.6$now,2)) > 0, round(bucket.2.6$now,2), NA), ifelse(length(round(rat.city$then.2.6,2)) > 0, round(rat.city$then.2.6,2), NA), ifelse(length(bucket.4.5$Actual.City) > 0, bucket.4.5$Actual.City, NA), ifelse(length(bucket.4.5$country) > 0, bucket.4.5$country, NA), ifelse(length(bucket.4.5$LATITUDE) > 0, bucket.4.5$LATITUDE, NA), ifelse(length(bucket.4.5$LONGITUDE) > 0, bucket.4.5$LONGITUDE, NA), ifelse(length(round(bucket.4.5$now,2)) > 0, round(bucket.4.5$now,2), NA), ifelse(length(round(rat.city$then.4.5,2)) > 0, round(rat.city$then.4.5,2), NA), ifelse(length(bucket.6.0$Actual.City) > 0, bucket.6.0$Actual.City, NA), ifelse(length(bucket.6.0$country) > 0, bucket.6.0$country, NA), ifelse(length(bucket.6.0$LATITUDE) > 0, bucket.6.0$LATITUDE, NA), ifelse(length(bucket.6.0$LONGITUDE) > 0, bucket.6.0$LONGITUDE, NA), ifelse(length(round(bucket.6.0$now,2)) > 0, round(bucket.6.0$now,2), NA), ifelse(length(round(rat.city$then.6.0,2)) > 0, round(rat.city$then.6.0,2), NA), ifelse(length(bucket.8.5$Actual.City) > 0, bucket.8.5$Actual.City, NA), ifelse(length(bucket.8.5$country) > 0, bucket.8.5$country, NA), ifelse(length(bucket.8.5$LATITUDE) > 0, bucket.8.5$LATITUDE, NA), ifelse(length(bucket.8.5$LONGITUDE) > 0, bucket.8.5$LONGITUDE, NA), ifelse(length(round(bucket.8.5$now,2)) > 0, round(bucket.8.5$now,2), NA), ifelse(length(round(rat.city$then.8.5,2)) > 0, round(rat.city$then.8.5,2), NA) )

	straight.dest.match <- rbind(straight.dest.match, one.dest.row)
}

colnames(straight.dest.match) <- c( "City.Name", "Country", "Latitude", "Longitude", "Temperature", "Compared.2.6.City.Name", "Compared.2.6.Country", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.City.Name", "Compared.4.5.Country", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.City.Name", "Compared.6.0.Country", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.City.Name", "Compared.8.5.Country", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature" )

write.csv(straight.dest.match, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_straightDestMatch.csv", row.names=FALSE)



##########################################################


for (i in seq(17.5, 44, by=0.5))
{
	ppr <- min( which ( bucket.points$lower > i))
	lwr <- max( which ( bucket.points$lower < i))
	print (paste(bucket.points$lower[lwr], "--" , i, "--" , bucket.points$lower[ppr]))
}


# bucket.cities <- merge(bucket.cities, country.codes, by.x="countryCode", by.y="code")
write.csv(bucket.cities, "~/Desktop/Tables/TheUpshot_ShiftingGlobalCities_buckets.csv", row.names=FALSE)


# TEST 20th CENTURY REANALYSIS GRID
rean20c.cities <- shift.cities
rean20c.JJA <- raster("/Users/dsmith/Desktop/Tables/shift cities/yycompos.kRxRqVQcWo.nc")
rean20c.DJF <- raster("/Users/dsmith/Desktop/Tables/shift cities/yycompos.eHH7o9Wjk2.nc")

rean20c.cities$rean20c[rean20c.cities$summerSzn=='JJA'] <- extract(rean20c.JJA,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear') - 273.15
rean20c.cities$rean20c[rean20c.cities$summerSzn=='DJF'] <- extract(rean20c.DJF,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear') - 273.15

rean20c.cities$model.2.6[rean20c.cities$summerSzn=='JJA'] <- extract(r.now.2.6.JJA,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')
rean20c.cities$model.2.6[rean20c.cities$summerSzn=='DJF'] <- extract(r.now.2.6.DJF,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

rean20c.cities$model.4.5[rean20c.cities$summerSzn=='JJA'] <- extract(r.now.4.5.JJA,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')
rean20c.cities$model.4.5[rean20c.cities$summerSzn=='DJF'] <- extract(r.now.4.5.DJF,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

rean20c.cities$model.6.0[rean20c.cities$summerSzn=='JJA'] <- extract(r.now.6.0.JJA,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')
rean20c.cities$model.6.0[rean20c.cities$summerSzn=='DJF'] <- extract(r.now.6.0.DJF,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')

rean20c.cities$model.8.5[rean20c.cities$summerSzn=='JJA'] <- extract(r.now.8.5.JJA,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='JJA',c("elong","LATITUDE")]), method='bilinear')
rean20c.cities$model.8.5[rean20c.cities$summerSzn=='DJF'] <- extract(r.now.8.5.DJF,  SpatialPoints(rean20c.cities[rean20c.cities$summerSzn=='DJF',c("elong","LATITUDE")]), method='bilinear')





plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$rean20c)
plot(rean20c.cities$now - rean20c.cities$rean20c, rean20c.cities$LATITUDE)

plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$model.2.6)
plot(rean20c.cities$now - rean20c.cities$model.2.6, rean20c.cities$LATITUDE)

plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$model.4.5)
plot(rean20c.cities$now - rean20c.cities$model.4.5, rean20c.cities$LATITUDE)

plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$model.6.0)
plot(rean20c.cities$now - rean20c.cities$model.6.0, rean20c.cities$LATITUDE)

plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$model.8.5)
plot(rean20c.cities$now - rean20c.cities$model.8.5, rean20c.cities$LATITUDE)


plot(rean20c.cities$diff.8.5, rean20c.cities$LATITUDE)
plot(rean20c.cities$now - rean20c.cities$model.8.5, rean20c.cities$LATITUDE)

plot(rean20c.cities$LONGITUDE, rean20c.cities$diff.8.5, ylim=c(-10.6700, 7.50))
plot(rean20c.cities$LONGITUDE, rean20c.cities$now - rean20c.cities$model.8.5, ylim=c(-10.6700, 7.50))

summary(rean20c.cities$now - rean20c.cities$model.8.5)
