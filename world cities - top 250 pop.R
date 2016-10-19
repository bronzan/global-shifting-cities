cities <- read.csv("~/Desktop/Tables/worldcitiespop.txt.gz")
cities.wPop <- subset(cities, !is.na(Population))
cities.top250 <- cities.wPop[order(cities.wPop$Population,decreasing=TRUE),][1:250,]
write.csv(cities.top250, "~/Desktop/Tables/top250worldCities.csv", row.names=FALSE)
cities.top101 <- cities.wPop[order(cities.wPop$Population,decreasing=TRUE),][1:101,]
write.csv(cities.top101, "~/Desktop/Tables/top101worldCities.csv", row.names=FALSE)

ffsta <- "//Users/dsmith/Desktop/Tables/ghcnd-stations.txt"
ffinv <- "//Users/dsmith/Desktop/Tables/ghcnd-inventory.txt"
fwsta <- c(11,-1,8,-1,9,-1,6,-1,2,-1,30,-1,3,-1,3,-1,5)
fwinv <- c(11,-1,8,-1,9,-1,4,-1,4,-1,4)
colssta <- c('ID','LATITUDE','LONGITUDE','ELEVATION','STATE','NAME','GSNFLAG','HCNFLAG','WMOID')
colsinv <- c('ID','LATITUDE','LONGITUDE','ELEMENT','FIRSTYEAR','LASTYEAR')
ALLsta <- read.fwf(ffsta, widths=fwsta, col.names=colssta, comment.char = "%")
# USsta <- ALLsta[which(substr(ALLsta$ID,1,2)=='US'),]
ALLinv <- read.fwf(ffinv, widths=fwinv, col.names=colsinv, comment.char = "%")
# USinv <- ALLinv[which(substr(ALLinv$ID,1,2)=='US'),]

library("geosphere")
library("curl")

ALLsta$close <- NA

for (i in 1:nrow(ALLsta))
{
  sta.ll <- ALLsta[i, c("LONGITUDE", "LATITUDE")]
  distances <- distHaversine(sta.ll, cities.top250[,c("Longitude","Latitude")])
  ALLsta$close[i] <- sum(distances <= 20000) > 0
}

closeStations <- subset(ALLsta, close==TRUE)
closeInventory <- subset(ALLinv, ELEMENT=='TMAX' & ID %in% closeStations$ID)
top250possibles <- merge(closeStations, closeInventory)
top250stations <- subset(top250possibles, FIRSTYEAR <= 1986 & LASTYEAR >= 2005)

write.csv(top250stations, "~/Desktop/Tables/top250stations.csv", row.names=FALSE)


fwdat <- c(11,4,2,4,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1)
colsdat <- c('ID', 'YEAR', 'MONTH', 'ELEMENT', 'VALUE1', 'MFLAG1', 'QFLAG1', 'SFLAG1', 'VALUE2', 'MFLAG2', 'QFLAG2', 'SFLAG2', 'VALUE3', 'MFLAG3', 'QFLAG3', 'SFLAG3', 'VALUE4', 'MFLAG4', 'QFLAG4', 'SFLAG4', 'VALUE5', 'MFLAG5', 'QFLAG5', 'SFLAG5', 'VALUE6', 'MFLAG6', 'QFLAG6', 'SFLAG6', 'VALUE7', 'MFLAG7', 'QFLAG7', 'SFLAG7', 'VALUE8', 'MFLAG8', 'QFLAG8', 'SFLAG8', 'VALUE9', 'MFLAG9', 'QFLAG9', 'SFLAG9', 'VALUE10', 'MFLAG10', 'QFLAG10', 'SFLAG10', 'VALUE11', 'MFLAG11', 'QFLAG11', 'SFLAG11', 'VALUE12', 'MFLAG12', 'QFLAG12', 'SFLAG12', 'VALUE13', 'MFLAG13', 'QFLAG13', 'SFLAG13', 'VALUE14', 'MFLAG14', 'QFLAG14', 'SFLAG14', 'VALUE15', 'MFLAG15', 'QFLAG15', 'SFLAG15', 'VALUE16', 'MFLAG16', 'QFLAG16', 'SFLAG16', 'VALUE17', 'MFLAG17', 'QFLAG17', 'SFLAG17', 'VALUE18', 'MFLAG18', 'QFLAG18', 'SFLAG18', 'VALUE19', 'MFLAG19', 'QFLAG19', 'SFLAG19', 'VALUE20', 'MFLAG20', 'QFLAG20', 'SFLAG20', 'VALUE21', 'MFLAG21', 'QFLAG21', 'SFLAG21', 'VALUE22', 'MFLAG22', 'QFLAG22', 'SFLAG22', 'VALUE23', 'MFLAG23', 'QFLAG23', 'SFLAG23', 'VALUE24', 'MFLAG24', 'QFLAG24', 'SFLAG24', 'VALUE25', 'MFLAG25', 'QFLAG25', 'SFLAG25', 'VALUE26', 'MFLAG26', 'QFLAG26', 'SFLAG26', 'VALUE27', 'MFLAG27', 'QFLAG27', 'SFLAG27', 'VALUE28', 'MFLAG28', 'QFLAG28', 'SFLAG28', 'VALUE29', 'MFLAG29', 'QFLAG29', 'SFLAG29', 'VALUE30', 'MFLAG30', 'QFLAG30', 'SFLAG30', 'VALUE31', 'MFLAG31', 'QFLAG31', 'SFLAG31')

top250stations$complete1 <- NA
top250stations$complete2 <- NA
stationData <- c()
for (j in 1:nrow(top250stations))
{
  print(paste0(j, " - ", date()))
  dataURL <- paste0("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/", top250stations$ID[j], ".dly")
  print(paste0(j, " - ", date(), dataURL))
  stationData <- c()
  rm(stationData)
  try(stationData <- read.fwf(curl(dataURL), widths=fwdat, col.names=colsdat, comment.char = "%", na.strings = "-9999"))
  if (top250stations$LATITUDE[j] <= 0) {
    stationTmax <- subset(stationData, ELEMENT=='TMAX' & YEAR %in% 1986:2005 & MONTH %in% c(12,1,2))
  } else {
    stationTmax <- subset(stationData, ELEMENT=='TMAX' & YEAR %in% 1986:2005 & MONTH %in% 6:8)
  }
  stationTmax$cover <- rowSums(is.na(stationTmax[,paste0("VALUE",1:28)]))
  stationTmax$complete <- stationTmax$cover/28 < 0.1
  top250stations$complete1[j] <- mean(is.na(stationTmax[,paste0("VALUE",1:28)])) < 0.1
  top250stations$complete2[j] <- mean(stationTmax$complete) > 0.9
}

write.csv(subset(top250stations, complete1==TRUE), "~/Desktop/Tables/top250cities-goodStations.csv", row.names=FALSE)
