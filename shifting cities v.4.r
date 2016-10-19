library(maps)
library(fields)
library(ncdf)
library(raster)
library(mapproj)
library(geosphere)
data(us.cities)

# cities.now.8.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_1986-2005/bcsd5/Extraction_tas_SummerAvgs.nc"
# cities.then.8.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_2081-2099/bcsd5/Extraction_tas_SummerAvgs.nc"
#
# cities.now.8.5.check.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_1986-2005_TAS/bcsd5/Extraction_tas_SummerAvgs.nc"
# cities.then.8.5.check.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_2081-2099_TAS/bcsd5/Extraction_tas_SummerAvgs.nc"
#
# # Now open the file and read its data
# cmip5.now.8.5.con <- open.ncdf(cities.now.8.5.file)
# cmip5.now.8.5.data <- get.var.ncdf(cmip5.now.8.5.con)
# close.ncdf(cmip5.now.8.5.con)
# cmip5.now.8.5.check.con <- open.ncdf(cities.now.8.5.check.file)
# cmip5.now.8.5.check.data <- get.var.ncdf(cmip5.now.8.5.check.con)
# close.ncdf(cmip5.now.8.5.check.con)
#
# cmip5.then.8.5.con <- open.ncdf(cities.then.8.5.file)
# cmip5.then.8.5.data <- get.var.ncdf(cmip5.then.8.5.con)
# close.ncdf(cmip5.then.8.5.con)
# cmip5.then.8.5.check.con <- open.ncdf(cities.then.8.5.check.file)
# cmip5.then.8.5.check.data <- get.var.ncdf(cmip5.then.8.5.check.con)
# close.ncdf(cmip5.then.8.5.check.con)
#
# cmip5.diff.8.5 <- cmip5.then.8.5.data - cmip5.now.8.5.data
# cmip5.diff.check.8.5 <- cmip5.then.8.5.check.data - cmip5.now.8.5.check.data
#
# #check against "old" TAS version
# mean(cmip5.diff.8.5 - cmip5.diff.check.8.5, na.rm=T)
#
#
#
# image.plot(cmip5.diff.8.5, zlim=c(min(cmip5.diff.8.5, na.rm=T), max(cmip5.diff.8.5, na.rm=T)), asp=.65, axes = FALSE, legend.shrink=0.5)

# START OVER USING PRISM AND RASTER PACKAGE
# Everything is in degrees Celsius, convert at the end

# FROM CMIP5
cities.now.2.6.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP26_allmonths_1986-2005/bcsd5/Extraction_tas_SummerAvgs.nc"
cities.then.2.6.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP26_allmonths_2081-2099/bcsd5/Extraction_tas_SummerAvgs.nc"

cities.now.4.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP45_allmonths_1986-2005/bcsd5/Extraction_tas_SummerAvgs.nc"
cities.then.4.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP45_allmonths_2081-2099/bcsd5/Extraction_tas_SummerAvgs.nc"

cities.now.6.0.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP60_allmonths_1986-2005/bcsd5/Extraction_tas_SummerAvgs.nc"
cities.then.6.0.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP60_allmonths_2081-2099/bcsd5/Extraction_tas_SummerAvgs.nc"

cities.now.8.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_1986-2005/bcsd5/Extraction_tas_SummerAvgs.nc"
cities.then.8.5.file <- "/Volumes/CCDatStat/CMIP5/CMIP5_RCP85_allmonths_2081-2099/bcsd5/Extraction_tas_SummerAvgs.nc"

r.now.2.6 <- raster(cities.now.2.6.file)
r.then.2.6 <- raster(cities.then.2.6.file)
r.diff.2.6 <- r.then.2.6-r.now.2.6

r.now.4.5 <- raster(cities.now.4.5.file)
r.then.4.5 <- raster(cities.then.4.5.file)
r.diff.4.5 <- r.then.4.5-r.now.4.5

r.now.6.0 <- raster(cities.now.6.0.file)
r.then.6.0 <- raster(cities.then.6.0.file)
r.diff.6.0 <- r.then.6.0-r.now.6.0

r.now.8.5 <- raster(cities.now.8.5.file)
r.then.8.5 <- raster(cities.then.8.5.file)
r.diff.8.5 <- r.then.8.5-r.now.8.5

# r.diff.8.52 <- raster("/Users/dsmith/Desktop/CMIP5_RCP85_tas_SummerAvgs_diff.nc")

# FROM PRISM
prism.path <- '/Volumes/CCDatStat/prism/monthly/tmax/1986/PRISM_tmax_stable_4kmM2_198606_bil/PRISM_tmax_stable_4kmM2_198606_bil.bil'
r.prism <- raster(prism.path)
r.prism <- r.prism-r.prism
for (yr in 1986:2005){
	for (mth in c("06","07","08")){
		prism.path <- paste0('/Volumes/CCDatStat/prism/monthly/tmax/', yr, '/PRISM_tmax_stable_4kmM2_', yr, mth, '_bil/PRISM_tmax_stable_4kmM2_', yr, mth, '_bil.bil')
		# add all 60 months together
		r.prism <- r.prism + raster(prism.path)
	}
}
# Divide by 60 to get average summer temp for 1986-2005
r.prism <- r.prism/60


# check all "markets" are in city list
tvm.sta <- read.csv("~/Desktop/Tables/Climate Matters weather station list .csv")
tvm.loc <- read.csv("~/Desktop/Tables/TV Mets locations 2014-07-01.csv")
tvm.loc$name <-  paste(tvm.loc$City, tvm.loc$State)
tvm.loc$Latitude <- round(tvm.loc$Latitude,2)
tvm.loc$Longitude <- round(tvm.loc$Longitude,2)

tvm.test <- merge(tvm.loc, us.cities, by="name", all.x=TRUE)
subset(tvm.test, is.na(pop))

# look at census.gov for a different list of cities...
# http://www.census.gov/popest/data/cities/totals/2013/SUB-EST2013.html
# http://www.census.gov/popest/data/cities/totals/2013/files/SUB-EST2013_ALL.csv
# cg.cities <- read.csv("~/Desktop/Tables/shift cities/census_gov pop SUB-EST2013_ALL.csv")
cg.cities <- read.csv(file = "http://www.census.gov/popest/data/cities/totals/2013/files/SUB-EST2013_ALL.csv")

dim(subset(cg.cities, CENSUS2010POP >= 40000 & SUMLEV %in% c(162, 170, 172)))
dim(subset(cg.cities, POPESTIMATE2010 >= 40000 & SUMLEV %in% c(162, 170, 172)))
dim(subset(cg.cities, POPESTIMATE2011 >= 40000 & SUMLEV %in% c(162, 170, 172)))
dim(subset(cg.cities, POPESTIMATE2012 >= 40000 & SUMLEV %in% c(162, 170, 172)))
dim(subset(cg.cities, POPESTIMATE2013 >= 40000 & SUMLEV %in% c(162, 170, 172)))

dim(subset(cg.cities, POPESTIMATE2010 >= 40000 & FUNCSTAT=="A"))
dim(subset(cg.cities, POPESTIMATE2011 >= 40000 & FUNCSTAT=="A"))
dim(subset(cg.cities, POPESTIMATE2012 >= 40000 & FUNCSTAT=="A"))
dim(subset(cg.cities, POPESTIMATE2013 >= 40000 & FUNCSTAT=="A"))

dim(subset(cg.cities, POPESTIMATE2013 >= 40000 & SUMLEV %in% c(162, 170, 172) & FUNCSTAT=="A"))

# Remove from destination list...
#	Lake Elsinore, CA
#	Pasadena, CA


# Change San Buenaventura to Ventura

# Send Utica to Texas?



shift.cities <- subset(us.cities, !(country.etc %in% c("HI", "AK")))
shift.cities$name <- substr(shift.cities$name, 1, nchar(shift.cities$name)-3)
shift.cities$elong <- 360 + shift.cities$long
tvm.loc$eLongitude <- 360 + tvm.loc$Longitude

# this is for handling NA projections for 10 southern cities for interactive
shift.cities$lat.plus1 <- shift.cities$lat + 0.125
shift.cities$lat.plus2 <- shift.cities$lat + 0.25
tvm.loc$Latitude.plus1 <- tvm.loc$Latitude + 0.125
tvm.loc$Latitude.plus2 <- tvm.loc$Latitude + 0.25

# 1.  find temperatures in PRISM data (r)
shift.cities$now <- extract(r.prism,  SpatialPoints(shift.cities[,c("long","lat")]), method='bilinear')
tvm.loc$now <- extract(r.prism,  SpatialPoints(tvm.loc[,c("Longitude","Latitude")]), method='bilinear')

# 2.  find model anomalies in CMIP5 data (r.diff)
shift.cities$diff.2.6 <- extract(r.diff.2.6,  SpatialPoints(shift.cities[,c("elong","lat")]), method='bilinear')
tvm.loc$diff.2.6 <- extract(r.diff.2.6,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude")]), method='bilinear')

shift.cities$diff.4.5 <- extract(r.diff.4.5,  SpatialPoints(shift.cities[,c("elong","lat")]), method='bilinear')
tvm.loc$diff.4.5 <- extract(r.diff.4.5,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude")]), method='bilinear')

shift.cities$diff.6.0 <- extract(r.diff.6.0,  SpatialPoints(shift.cities[,c("elong","lat")]), method='bilinear')
tvm.loc$diff.6.0 <- extract(r.diff.6.0,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude")]), method='bilinear')

shift.cities$diff.8.5 <- extract(r.diff.8.5,  SpatialPoints(shift.cities[,c("elong","lat")]), method='bilinear')
# more NA handling
diff.8.5.plus1 <- extract(r.diff.8.5,  SpatialPoints(shift.cities[,c("elong","lat.plus1")]), method='bilinear')
diff.8.5.plus2 <- extract(r.diff.8.5,  SpatialPoints(shift.cities[,c("elong","lat.plus2")]), method='bilinear')
shift.cities$diff.8.5[is.na(shift.cities$diff.8.5)] <- diff.8.5.plus1[is.na(shift.cities$diff.8.5)]
shift.cities$diff.8.5[is.na(shift.cities$diff.8.5)] <- diff.8.5.plus2[is.na(shift.cities$diff.8.5)]

tvm.loc$diff.8.5 <- extract(r.diff.8.5,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude")]), method='bilinear')
# more NA handling
tvm.diff.8.5.plus1 <- extract(r.diff.8.5,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude.plus1")]), method='bilinear')
tvm.diff.8.5.plus2 <- extract(r.diff.8.5,  SpatialPoints(tvm.loc[,c("eLongitude","Latitude.plus2")]), method='bilinear')
tvm.loc$diff.8.5[is.na(tvm.loc$diff.8.5)] <- tvm.diff.8.5.plus1[is.na(tvm.loc$diff.8.5)]
tvm.loc$diff.8.5[is.na(tvm.loc$diff.8.5)] <- tvm.diff.8.5.plus2[is.na(tvm.loc$diff.8.5)]

# 3.  add anomolies to temps
shift.cities$then.2.6 <- shift.cities$now + shift.cities$diff.2.6
tvm.loc$then.2.6 <- tvm.loc$now + tvm.loc$diff.2.6
shift.cities$then.4.5 <- shift.cities$now + shift.cities$diff.4.5
tvm.loc$then.4.5 <- tvm.loc$now + tvm.loc$diff.4.5
shift.cities$then.6.0 <- shift.cities$now + shift.cities$diff.6.0
tvm.loc$then.6.0 <- tvm.loc$now + tvm.loc$diff.6.0
shift.cities$then.8.5 <- shift.cities$now + shift.cities$diff.8.5
tvm.loc$then.8.5 <- tvm.loc$now + tvm.loc$diff.8.5

us.capitals <- subset(shift.cities, capital != 0)

# Luxor, Egypt: 41.1 (on par w/ Lake Havasu, AZ)
# Khartoum, Sudan: 41.7
# Abu Dhabi, UAE: 43
# Riyadh, Saudi Arabia: 43.9
# And...Kuwait City, Kuwait: 45.6

mideast <- data.frame(rbind(c("Luxor", "Egypt", 25.70, 32.65,41.1),
				c("Khartoum", "Sudan", 15.58, 32.52, 41.7),
				c("Abu Dhabi", "UAE", 24.48, 54.37, 43),
				c("Riyadh", "Saudi Arabia", 24.65, 46.77, 43.9),
				c("Kuwait City", "Kuwait", 29.38, 47.99, 45.6)))
colnames(mideast) <- c("name", "country.etc", "lat", "long", "now")
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


match.cities <- subset(shift.cities, !(name %in% c("Lake Elsinore", "Pasadena") & country.etc == "CA") )

all.city.matches.out <- data.frame()
radish.matches.out <- data.frame()
all.matches.out <- data.frame()
# 4.  capitals loop
# for (i in 1:nrow(us.capitals)){
for (i in 1:nrow(shift.cities)){

	# capital <- us.capitals[i,]
	capital <- shift.cities[i,]

# 	a.  find other cities that are SSE-S-SSW of current city
	cap.bearing <- bearing(capital[,c("long","lat")], match.cities[,c("long","lat")])

# 	b.  calculate distances
	cap.distance <- distHaversine(capital[,c("long","lat")], match.cities[,c("long","lat")], r=6378137)

#	b.2. calculate latitude and longitude differences
	cap.lat.diff <-  capital$lat - match.cities$lat
	cap.long.diff <- capital$long - match.cities$long

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

	temp.order.long.diff.2.6 <- cap.long.diff[order(abs(cap.diff.2.6), decreasing=FALSE)]
	temp.order.long.diff.4.5 <- cap.long.diff[order(abs(cap.diff.4.5), decreasing=FALSE)]
	temp.order.long.diff.6.0 <- cap.long.diff[order(abs(cap.diff.6.0), decreasing=FALSE)]
	temp.order.long.diff.8.5 <- cap.long.diff[order(abs(cap.diff.8.5), decreasing=FALSE)]

	# closest.temp.bearing.2.6 <- abs(temp.order.bearing.2.6[1:10]-180)
	# closest.temp.bearing.4.5 <- abs(temp.order.bearing.4.5[1:10]-180)
	# closest.temp.bearing.6.0 <- abs(temp.order.bearing.6.0[1:10]-180)
	# closest.temp.bearing.8.5 <- abs(temp.order.bearing.8.5[1:10]-180)

	closest.temp.bearing.2.6 <- abs( abs(temp.order.bearing.2.6[1:10]) - 180 )
	closest.temp.bearing.4.5 <- abs( abs(temp.order.bearing.4.5[1:10]) - 180 )
	closest.temp.bearing.6.0 <- abs( abs(temp.order.bearing.6.0[1:10]) - 180 )
	closest.temp.bearing.8.5 <- abs( abs(temp.order.bearing.8.5[1:10]) - 180 )

### NEW : TOP 10 CLOSEST BY TEMP
###		  TOP 3 BEST BEARING (closest to southern)
###		  TOP 1 MOST SOUTH (biggest lat diff)

	bearing.three.2.6 <- temp.order.cities.2.6[1:10,][order(closest.temp.bearing.2.6, decreasing=FALSE),][1:3,]
	bearing.three.4.5 <- temp.order.cities.4.5[1:10,][order(closest.temp.bearing.4.5, decreasing=FALSE),][1:3,]
	bearing.three.6.0 <- temp.order.cities.6.0[1:10,][order(closest.temp.bearing.6.0, decreasing=FALSE),][1:3,]
	bearing.three.8.5 <- temp.order.cities.8.5[1:10,][order(closest.temp.bearing.8.5, decreasing=FALSE),][1:3,]

	lat.diff.three.2.6 <- temp.order.lat.diff.2.6[1:10][order(closest.temp.bearing.2.6, decreasing=FALSE)][1:3]
	lat.diff.three.4.5 <- temp.order.lat.diff.4.5[1:10][order(closest.temp.bearing.4.5, decreasing=FALSE)][1:3]
	lat.diff.three.6.0 <- temp.order.lat.diff.6.0[1:10][order(closest.temp.bearing.6.0, decreasing=FALSE)][1:3]
	lat.diff.three.8.5 <- temp.order.lat.diff.8.5[1:10][order(closest.temp.bearing.8.5, decreasing=FALSE)][1:3]

	mtch.2.6 <- bearing.three.2.6[which.max(lat.diff.three.2.6),]
	mtch.4.5 <- bearing.three.4.5[which.max(lat.diff.three.4.5),]
	mtch.6.0 <- bearing.three.6.0[which.max(lat.diff.three.6.0),]
	mtch.8.5 <- bearing.three.8.5[which.max(lat.diff.three.8.5),]


	print(paste0(capital$name, ", ", capital$country.etc, " - ", mtch.2.6$name, ", ", mtch.2.6$country.etc, " - ", mtch.4.5$name, ", ", mtch.4.5$country.etc, " - ", mtch.6.0$name, ", ", mtch.6.0$country.etc, " - ", mtch.8.5$name, ", ", mtch.8.5$country.etc))

	# print(paste0(capital$name, ", ", capital$country.etc, " - ", dim(ordered.cities.2.6)[1], " - ", dim(ordered.cities.4.5)[1], " - ", dim(ordered.cities.6.0)[1], dim(ordered.cities.8.5)[1]))

	all.city.matches.out <- rbind(all.city.matches.out, c(city=capital$name, mtch.2.6=mtch.2.6$name, mtch.4.5=mtch.4.5$name, mtch.6.0=mtch.6.0$name, mtch.8.5=mtch.8.5$name))

	one.radish.row <- c(capital$name, capital$country.etc, capital$lat, capital$long, round(capital$now,2), mtch.8.5$name, mtch.8.5$country.etc, mtch.8.5$lat, mtch.8.5$long, round(mtch.8.5$now,2), round(capital$then.8.5,2))

	one.selcaps.row <- c(capital$name, capital$country.etc, capital$lat, capital$long, round(capital$now,2), mtch.2.6$name, mtch.2.6$country.etc, mtch.2.6$lat, mtch.2.6$long, round(mtch.2.6$now,2), round(capital$then.2.6,2), mtch.4.5$name, mtch.4.5$country.etc, mtch.4.5$lat, mtch.4.5$long, round(mtch.4.5$now,2), round(capital$then.4.5,2), mtch.6.0$name, mtch.6.0$country.etc, mtch.6.0$lat, mtch.6.0$long, round(mtch.6.0$now,2), round(capital$then.6.0,2), mtch.8.5$name, mtch.8.5$country.etc, mtch.8.5$lat, mtch.8.5$long, round(mtch.8.5$now,2), round(capital$then.8.5,2))

	if (!is.na(capital$then.8.5) & capital$then.8.5 > 41.48413) {
		me.diff <- abs(capital$then.8.5 - mideast$now)
		me.match <- mideast[which.min(me.diff),]

		one.radish.row <- c(capital$name, capital$country.etc, capital$lat, capital$long, round(capital$now,2), me.match$name, me.match$country.etc, me.match$lat, me.match$long, me.match$now, round(capital$then.8.5,2))

		one.selcaps.row <- c(capital$name, capital$country.etc, capital$lat, capital$long, round(capital$now,2), mtch.2.6$name, mtch.2.6$country.etc, mtch.2.6$lat, mtch.2.6$long, round(mtch.2.6$now,2), round(capital$then.2.6,2), mtch.4.5$name, mtch.4.5$country.etc, mtch.4.5$lat, mtch.4.5$long, round(mtch.4.5$now,2), round(capital$then.4.5,2), mtch.6.0$name, mtch.6.0$country.etc, mtch.6.0$lat, mtch.6.0$long, round(mtch.6.0$now,2), round(capital$then.6.0,2), me.match$name, me.match$country.etc, me.match$lat, me.match$long, round(me.match$now,2), round(capital$then.8.5,2))
	}
	radish.matches.out <- rbind(radish.matches.out, one.radish.row)
	all.matches.out <- rbind(all.matches.out, one.selcaps.row)
}

colnames(radish.matches.out) <- c("City.Name", "State", "Latitude", "Longitude", "Temperature", "Compared.City.Name", "Compared.State", "Compared.Latitude", "Compared.Longitude", "Compared.Temperature", "Projected.Temperature")
colnames(all.matches.out) <- c("City.Name", "State", "Latitude", "Longitude", "Temperature", "Compared.2.6.City.Name", "Compared.2.6.State", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.City.Name", "Compared.4.5.State", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.City.Name", "Compared.6.0.State", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.City.Name", "Compared.8.5.State", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature")

#write.csv(radish.matches.out, "~/Desktop/Tables/Radish_ShiftingCities_3stepMatch.csv", row.names=FALSE)
#write.csv(all.matches.out, "~/Desktop/Tables/TheUpshot_ShiftingCities_3stepMatch.csv", row.names=FALSE)

# check inconsistent 8.5 matches 2015-12-02
radish.matches.check <- read.csv("~/Desktop/Tables/shift cities/Radish_ShiftingCities_3stepMatch_final.csv")
radish.matches.check2 <- read.csv("~/Desktop/Tables/shift cities/Radish_ShiftingCities_3stepMatch_noNAs.csv")

all.matches.ventura <- all.matches.out
all.matches.ventura$City.Name[all.matches.ventura$City.Name == 'San Buenaventura'] <- 'Ventura'
all.graft.old <- merge(all.matches.ventura, radish.matches.check, by=c("City.Name","State"), all=TRUE)

# sum( all.graft.old[,c("Compared.8.5.City.Name")] == all.graft.old[,c("Compared.City.Name")] )
# sum( all.graft.old[,c("Compared.8.5.State")] == all.graft.old[,c("Compared.State")] )
# sum( all.graft.old[,c("Compared.8.5.Latitude")] == all.graft.old[,c("Compared.Latitude")] )
# sum( all.graft.old[,c("Compared.8.5.Longitude")] == all.graft.old[,c("Compared.Longitude")] )
# sum( all.graft.old[,c("Compared.8.5.Temperature")] == (all.graft.old[,c("Compared.Temperature")]-32)/1.8 )
# sum( all.graft.old[,c("Projected.8.5.Temperature")] == all.graft.old[,c("Projected.Temperature")] )

all.graft.old[,c("Compared.8.5.City.Name")] <- all.graft.old[,c("Compared.City.Name")]
all.graft.old[,c("Compared.8.5.State")] <- all.graft.old[,c("Compared.State")]
all.graft.old[,c("Compared.8.5.Latitude")] <- all.graft.old[,c("Compared.Latitude")]
all.graft.old[,c("Compared.8.5.Longitude")] <- all.graft.old[,c("Compared.Longitude")]
all.graft.old[,c("Compared.8.5.Temperature")] <- (all.graft.old[,c("Compared.Temperature")]-32)/1.8
all.graft.old[,c("Projected.8.5.Temperature")] <- all.graft.old[,c("Projected.Temperature")]
all.graft.old <- all.graft.old[,c("City.Name", "State", "Latitude.x", "Longitude.x", "Temperature.x", "Compared.2.6.City.Name", "Compared.2.6.State", "Compared.2.6.Latitude", "Compared.2.6.Longitude", "Compared.2.6.Temperature", "Projected.2.6.Temperature", "Compared.4.5.City.Name", "Compared.4.5.State", "Compared.4.5.Latitude", "Compared.4.5.Longitude", "Compared.4.5.Temperature", "Projected.4.5.Temperature", "Compared.6.0.City.Name", "Compared.6.0.State", "Compared.6.0.Latitude", "Compared.6.0.Longitude", "Compared.6.0.Temperature", "Projected.6.0.Temperature", "Compared.8.5.City.Name", "Compared.8.5.State", "Compared.8.5.Latitude", "Compared.8.5.Longitude", "Compared.8.5.Temperature", "Projected.8.5.Temperature")]
colnames(all.graft.old) <- c("City_Name", "State", "Latitude", "Longitude", "Temperature", "Compared_26_City_Name", "Compared_26_State", "Compared_26_Latitude", "Compared_26_Longitude", "Compared_26_Temperature", "Projected_26_Temperature", "Compared_45_City_Name", "Compared_45_State", "Compared_45_Latitude", "Compared_45_Longitude", "Compared_45_Temperature", "Projected_45_Temperature", "Compared_60_City_Name", "Compared_60_State", "Compared_60_Latitude", "Compared_60_Longitude", "Compared_60_Temperature", "Projected_60_Temperature", "Compared_85_City_Name", "Compared_85_State", "Compared_85_Latitude", "Compared_85_Longitude", "Compared_85_Temperature", "Projected_85_Temperature")

write.csv(all.graft.old, "~/Desktop/Tables/TheUpshot_ShiftingCities_3stepMatch.csv", row.names=FALSE)

head(radish.matches.out[ ,c("City.Name", "Compared.City.Name", "Projected.Temperature")])
head(radish.matches.check[ ,c("City.Name", "Compared.City.Name", "Projected.Temperature")])
head(radish.matches.check2[ ,c("City.Name", "Compared.City.Name", "Projected.Temperature")])

head(radish.matches.out[ ,c("City.Name", "Compared.City.Name", "Compared.Temperature")])
head(radish.matches.check[ ,c("City.Name", "Compared.City.Name", "Compared.Temperature")])
head(radish.matches.check2[ ,c("City.Name", "Compared.City.Name", "Compared.Temperature")])
head(radish.matches.check2$Compared.Temperature*1.8 + 32)

(radish.matches.out$Compared.City.Name==radish.matches.check$Compared.City.Name)[1:20]

sum(radish.matches.out$Compared.City.Name==radish.matches.check2$Compared.City.Name)






## OTHER MIDDLE EAST EXPLORATION CODE

me <- read.csv("~/Desktop/Tables/shift cities/shiftingcities-middleeast.csv")
me$TMAX[which(me$Quality.Flag %in% c("D","G","I","O"))] <- NA
me.means <- aggregate(TMAX ~ STATION + STATION_NAME, data=me, mean, na.rm=TRUE)
me.means$TMAX <- me.means$TMAX/10
me.counts <- data.frame(table(me$STATION))
colnames(me.counts) <- c("STATION", "FREQ")

middle.east <- merge(me.means, me.counts)
middle.east.out <- subset(middle.east, FREQ >= 1656)

write.csv(middle.east.out, "~/Desktop/Tables/shift cities/other possible middle east locations.csv")

#	limit to 101 buckets?
#		use min and max projected temp to filter current cities and then come up with buckets
# BUCKET CREATION
	cutoff <- min(shift.cities$then.8.5, na.rm=TRUE) - 0.5
	match.cities <- subset(shift.cities, now >= cutoff)
	match.cities$buckets <-  as.numeric(cut(match.cities$now,101))
	write.csv(match.cities, "~/Desktop/tables/shift cities/match buckets.csv", row.names=FALSE)

	pdf('~/Desktop/Plots/Shifting Cities Population filters.pdf',  width=11, height=7)
		map('state', col="white", fill = TRUE, resolution = 0, lty = 0, bg='light gray',)
		map('state', add=TRUE, names=TRUE)
		points(shift.cities$long, shift.cities$lat, col='black', bg='blue', pch=21, cex=.7)
		title("1001 cities > 40k pop", col.main = "black")

		map('state', col="white", fill = TRUE, resolution = 0, lty = 0, bg='light gray')
		map('state', add=TRUE, names=TRUE)
		points(subset(match.cities, pop > 100000)$long, subset(match.cities, pop > 100000)$lat, col='black', bg='green', pch=21, cex=.7)
		title("258 match cities > 100k pop", col.main = "black")

		map('state', col="white", fill = TRUE, resolution = 0, lty = 0, bg='light gray')
		map('state', add=TRUE, names=TRUE)
		points(subset(match.cities, pop > 200000)$long, subset(match.cities, pop > 200000)$lat, col='black', bg='red', pch=21, cex=.7)
		title("94 match cities > 200k pop", col.main = "black")
	dev.off()




tvmets.matches.out <- data.frame()
# 4.  tv mets loop
for (i in 1:nrow(tvm.loc)){
	market <- tvm.loc[i,]

# 	a.  find other cities that are SSE-S-SSW of current city
	mkt.bearing <- bearing(market[,c("Longitude","Latitude")], match.cities[,c("long","lat")])

# 	b.  calculate distances
	mkt.distance <- distHaversine(market[,c("Longitude","Latitude")], match.cities[,c("long","lat")], r=6378137)

#	b.2. calculate Latitude and Longitude differences
	mkt.Latitude.diff <-  market$Latitude - match.cities$lat
	mkt.Longitude.diff <- market$Longitude - match.cities$long

# 	c.  calculate temperature difference
	mkt.diff.2.6 <- market$then.2.6 - match.cities$now
	mkt.diff.4.5 <- market$then.4.5 - match.cities$now
	mkt.diff.6.0 <- market$then.6.0 - match.cities$now
	mkt.diff.8.5 <- market$then.8.5 - match.cities$now

#	d.  select all with temp match within 0.5 degrees then look at "best" bearing

	temp.order.cities.2.6 <- match.cities[order(abs(mkt.diff.2.6), decreasing=FALSE),]
	temp.order.cities.4.5 <- match.cities[order(abs(mkt.diff.4.5), decreasing=FALSE),]
	temp.order.cities.6.0 <- match.cities[order(abs(mkt.diff.6.0), decreasing=FALSE),]
	temp.order.cities.8.5 <- match.cities[order(abs(mkt.diff.8.5), decreasing=FALSE),]

	temp.order.bearing.2.6 <- mkt.bearing[order(abs(mkt.diff.2.6), decreasing=FALSE)]
	temp.order.bearing.4.5 <- mkt.bearing[order(abs(mkt.diff.4.5), decreasing=FALSE)]
	temp.order.bearing.6.0 <- mkt.bearing[order(abs(mkt.diff.6.0), decreasing=FALSE)]
	temp.order.bearing.8.5 <- mkt.bearing[order(abs(mkt.diff.8.5), decreasing=FALSE)]

	temp.order.Latitude.diff.2.6 <- mkt.Latitude.diff[order(abs(mkt.diff.2.6), decreasing=FALSE)]
	temp.order.Latitude.diff.4.5 <- mkt.Latitude.diff[order(abs(mkt.diff.4.5), decreasing=FALSE)]
	temp.order.Latitude.diff.6.0 <- mkt.Latitude.diff[order(abs(mkt.diff.6.0), decreasing=FALSE)]
	temp.order.Latitude.diff.8.5 <- mkt.Latitude.diff[order(abs(mkt.diff.8.5), decreasing=FALSE)]

	temp.order.Longitude.diff.2.6 <- mkt.Longitude.diff[order(abs(mkt.diff.2.6), decreasing=FALSE)]
	temp.order.Longitude.diff.4.5 <- mkt.Longitude.diff[order(abs(mkt.diff.4.5), decreasing=FALSE)]
	temp.order.Longitude.diff.6.0 <- mkt.Longitude.diff[order(abs(mkt.diff.6.0), decreasing=FALSE)]
	temp.order.Longitude.diff.8.5 <- mkt.Longitude.diff[order(abs(mkt.diff.8.5), decreasing=FALSE)]

	closest.temp.bearing.2.6 <- abs(temp.order.bearing.2.6[1:10]-180)
	closest.temp.bearing.4.5 <- abs(temp.order.bearing.4.5[1:10]-180)
	closest.temp.bearing.6.0 <- abs(temp.order.bearing.6.0[1:10]-180)
	closest.temp.bearing.8.5 <- abs(temp.order.bearing.8.5[1:10]-180)


### NEW : TOP 10 CLOSEST BY TEMP
###		  TOP 3 BEST BEARING (closest to southern)
###		  TOP 1 MOST SOUTH (biggest Latitude diff)

	bearing.three.2.6 <- temp.order.cities.2.6[1:10,][order(closest.temp.bearing.2.6, decreasing=FALSE),][1:3,]
	bearing.three.4.5 <- temp.order.cities.4.5[1:10,][order(closest.temp.bearing.4.5, decreasing=FALSE),][1:3,]
	bearing.three.6.0 <- temp.order.cities.6.0[1:10,][order(closest.temp.bearing.6.0, decreasing=FALSE),][1:3,]
	bearing.three.8.5 <- temp.order.cities.8.5[1:10,][order(closest.temp.bearing.8.5, decreasing=FALSE),][1:3,]

	Latitude.diff.three.2.6 <- temp.order.Latitude.diff.2.6[1:10][order(closest.temp.bearing.2.6, decreasing=FALSE)][1:3]
	Latitude.diff.three.4.5 <- temp.order.Latitude.diff.4.5[1:10][order(closest.temp.bearing.4.5, decreasing=FALSE)][1:3]
	Latitude.diff.three.6.0 <- temp.order.Latitude.diff.6.0[1:10][order(closest.temp.bearing.6.0, decreasing=FALSE)][1:3]
	Latitude.diff.three.8.5 <- temp.order.Latitude.diff.8.5[1:10][order(closest.temp.bearing.8.5, decreasing=FALSE)][1:3]

	mtch.2.6 <- bearing.three.2.6[which.max(Latitude.diff.three.2.6),]
	mtch.4.5 <- bearing.three.4.5[which.max(Latitude.diff.three.4.5),]
	mtch.6.0 <- bearing.three.6.0[which.max(Latitude.diff.three.6.0),]
	mtch.8.5 <- bearing.three.8.5[which.max(Latitude.diff.three.8.5),]


	print(paste0(market$Geographic.location, ", ", market$State, " - ", mtch.2.6$name, " - ", mtch.4.5$name, " - ", mtch.6.0$name, " - ", mtch.8.5$name))

	# print(paste0(market$name, ", ", market$country.etc, " - ", dim(ordered.cities.2.6)[1], " - ", dim(ordered.cities.4.5)[1], " - ", dim(ordered.cities.6.0)[1], dim(ordered.cities.8.5)[1]))

# 	all.city.matches.out <- rbind(all.city.matches.out, c(city=market$name, mtch.2.6=mtch.2.6$name, mtch.4.5=mtch.4.5$name, mtch.6.0=mtch.6.0$name, mtch.8.5=mtch.8.5$name))
	one.tvmets.row <- c(market$Geographic.location, market$State, market$Latitude, market$Longitude, round(market$now,2), mtch.8.5$name, mtch.8.5$country.etc, mtch.8.5$lat, mtch.8.5$long, round(mtch.8.5$now,2), round(market$then.8.5,2))
	if (!is.na(market$then.8.5) & market$then.8.5 > 41.48413) {
		me.diff <- abs(market$then.8.5 - mideast$now)
		me.match <- mideast[which.min(me.diff),]
		one.tvmets.row <- c(market$Geographic.location, market$State, market$Latitude, market$Longitude, round(market$now,2), me.match$name, me.match$country.etc, me.match$lat, me.match$long, me.match$now, round(market$then.8.5,2))

	}

	tvmets.matches.out <- rbind(tvmets.matches.out, one.tvmets.row)
}

colnames(tvmets.matches.out) <- c("City.Name", "State", "Latitude", "Longitude", "Temperature", "Compared.City.Name", "Compared.State", "Compared.Latitude", "Compared.Longitude", "Compared.Temperature", "Projected.Temperature")

 write.csv(tvmets.matches.out, "~/Desktop/Tables/TVMets_ShiftingCities_3stepMatch.csv", row.names=FALSE)
