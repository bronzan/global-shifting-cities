library("curl")

stationList <- read.csv("~/Desktop/Tables/GlobalCitiesList.csv")

fwdat <- c(11,4,2,4,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1,5,1,1,1)
colsdat <- c('ID', 'YEAR', 'MONTH', 'ELEMENT', 'VALUE1', 'MFLAG1', 'QFLAG1', 'SFLAG1', 'VALUE2', 'MFLAG2', 'QFLAG2', 'SFLAG2', 'VALUE3', 'MFLAG3', 'QFLAG3', 'SFLAG3', 'VALUE4', 'MFLAG4', 'QFLAG4', 'SFLAG4', 'VALUE5', 'MFLAG5', 'QFLAG5', 'SFLAG5', 'VALUE6', 'MFLAG6', 'QFLAG6', 'SFLAG6', 'VALUE7', 'MFLAG7', 'QFLAG7', 'SFLAG7', 'VALUE8', 'MFLAG8', 'QFLAG8', 'SFLAG8', 'VALUE9', 'MFLAG9', 'QFLAG9', 'SFLAG9', 'VALUE10', 'MFLAG10', 'QFLAG10', 'SFLAG10', 'VALUE11', 'MFLAG11', 'QFLAG11', 'SFLAG11', 'VALUE12', 'MFLAG12', 'QFLAG12', 'SFLAG12', 'VALUE13', 'MFLAG13', 'QFLAG13', 'SFLAG13', 'VALUE14', 'MFLAG14', 'QFLAG14', 'SFLAG14', 'VALUE15', 'MFLAG15', 'QFLAG15', 'SFLAG15', 'VALUE16', 'MFLAG16', 'QFLAG16', 'SFLAG16', 'VALUE17', 'MFLAG17', 'QFLAG17', 'SFLAG17', 'VALUE18', 'MFLAG18', 'QFLAG18', 'SFLAG18', 'VALUE19', 'MFLAG19', 'QFLAG19', 'SFLAG19', 'VALUE20', 'MFLAG20', 'QFLAG20', 'SFLAG20', 'VALUE21', 'MFLAG21', 'QFLAG21', 'SFLAG21', 'VALUE22', 'MFLAG22', 'QFLAG22', 'SFLAG22', 'VALUE23', 'MFLAG23', 'QFLAG23', 'SFLAG23', 'VALUE24', 'MFLAG24', 'QFLAG24', 'SFLAG24', 'VALUE25', 'MFLAG25', 'QFLAG25', 'SFLAG25', 'VALUE26', 'MFLAG26', 'QFLAG26', 'SFLAG26', 'VALUE27', 'MFLAG27', 'QFLAG27', 'SFLAG27', 'VALUE28', 'MFLAG28', 'QFLAG28', 'SFLAG28', 'VALUE29', 'MFLAG29', 'QFLAG29', 'SFLAG29', 'VALUE30', 'MFLAG30', 'QFLAG30', 'SFLAG30', 'VALUE31', 'MFLAG31', 'QFLAG31', 'SFLAG31')

stationList$summerMean1 <- NA
stationList$summerMean2 <- NA
stationData <- c()
for (j in 1:nrow(stationList))
{
  print(paste0(j, " - ", date()))
  dataURL <- paste0("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/", stationList$ID[j], ".dly")
  print(paste0(j, " - ", date(), dataURL))
  stationData <- c()
  rm(stationData)
  try()

  if (stationList$LATITUDE[j] <= 0) {
    stationTmax <- subset(stationData, ELEMENT=='TMAX' & YEAR %in% 1986:2005 & MONTH %in% c(12,1,2))
  } else {
    stationTmax <- subset(stationData, ELEMENT=='TMAX' & YEAR %in% 1986:2005 & MONTH %in% 6:8)
  }
  stationTmax$moMean <- rowMean(stationTmax[,paste0("VALUE",1:31)], na.rm=TRUE)
  stationList$summerMean1[j] <- mean(stationTmax$moMean)
  stationList$summerMean2[j] <- mean(stationTmax[,paste0("VALUE",1:31)], na.rm=TRUE)
}


# USING GHCN-M INSTEAD?
library("geosphere")
ffmax <- "~/Desktop/Tables/ghcnm.v3.3.0.20151128/ghcnm.tmax.v3.3.0.20151128.qca.dat"
ffinv <- "~/Desktop/Tables/ghcnm.v3.3.0.20151128/ghcnm.tmax.v3.3.0.20151128.qca.inv"
# fwdat<-c(6,1,4,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1,-1,5,1)
fwdat <- c(11, 4, 4, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1)
fwinv <- c(11, -1, 8, -1, 9, -1, 6, -1, 30, -1, 4, 1, 5, 2, 2, 2, 2, 1, 2, 16, 1)

coldat <- c('ID','YEAR','ELEMENT','Jan', 'dmJan', 'qcJan', 'dsJan', 'Feb', 'dmFeb', 'qcFeb', 'dsFeb', 'Mar', 'dmMar', 'qcMar', 'dsMar', 'Apr', 'dmApr', 'qcApr', 'dsApr', 'May', 'dmMay', 'qcMay', 'dsMay',  'Jun', 'dmJun', 'qcJun', 'dsJun', 'Jul', 'dmJul', 'qcJul', 'dsJul', 'Aug', 'dmAug', 'qcAug', 'dsAug', 'Sep', 'dmSep', 'qcSep', 'dsSep', 'Oct', 'dmOct', 'qcOct', 'dsOct', 'Nov', 'dmNov', 'qcNov', 'dsNov', 'Dec', 'dmDec', 'qcDec', 'dsDec')

colinv <- c("ID", "LATITUDE", "LONGITUDE", "STNELEV", "NAME", "GRELEV", "POPCLS", "POPSIZ", "TOPO", "STVEG", "STLOC", "OCNDIS", "AIRSTN", "TOWNDIS", "GRVEG", "POPCSS")

ghcnm.max.inv <- read.fwf(ffinv, widths=fwinv, col.names=colinv, comment.char = "%", na.strings="-9999")
ghcnm.max.inv$countryCode <- substr(ghcnm.max.inv$ID, 1, 3)
ghcnm.max.inv$WMOID <- substr(ghcnm.max.inv$ID, 4, 8)

ghcnm.max <- read.fwf(ffmax, widths=fwdat, col.names=coldat, comment.char = "%", na.strings="-9999")
ghcnm.max.t <- subset(ghcnm.max, YEAR >= 1986 & YEAR <= 2005)[,c("ID", "YEAR", "Jan", "Feb", "Jun", "Jul", "Aug", "Dec")]
ghcnm.max.t$countryCode <- substr(ghcnm.max.t$ID,1,3)
ghcnm.max.t$WMOID <- substr(ghcnm.max.t$ID,4,8)

cities <- read.csv("~/Desktop/Tables/worldcitiespop.txt.gz")
cities.wPop <- subset(cities, !is.na(Population))
cities.top250 <- cities.wPop[order(cities.wPop$Population,decreasing=TRUE),][1:250,]
cities.top101 <- cities.wPop[order(cities.wPop$Population,decreasing=TRUE),][1:101,]

for (i in 1:nrow(ghcnm.max.inv))
{
  sta.ll <- ghcnm.max.inv[i, c("LONGITUDE", "LATITUDE")]
  distances <- distHaversine(sta.ll, cities.top250[,c("Longitude","Latitude")])
  ghcnm.max.inv$close[i] <- sum(distances <= 20000) > 0
}

closeStations <- subset(ghcnm.max.inv, close==TRUE)

closeStations <- merge(closeStations, ghcnm.max.t, by="ID")

wmoidStations <- unique(closeStations$ID)
results <- c()
for (id in wmoidStations){
   one.station <- subset(closeStations, ID == id)
   if (nrow(one.station) >= 0.9){
      summer <- c("Jun", "Jul", "Aug")
      summerSzn <- "JJA"
      if (one.station$LATITUDE[1] <= 0){
         summer <- c("Dec", "Jan", "Feb")
         summerSzn <- "DJF"
      }
      one.summers <- one.station[,summer]
      summerMean <- round(mean(as.matrix(one.summers), na.rm=TRUE)/100,2)
      one.result <- cbind(one.station[1, c("ID", "LATITUDE", "LONGITUDE", "NAME")], summerMean, summerSzn)
      results <- rbind(results, one.result)
   }
}

results$countryCode <- substr(results$ID,1,3)
results$WMOID <- substr(results$ID,4,8)
write.csv(results, "~/Desktop/Tables/globalSummerMeans.csv", row.names=FALSE)




#file handling for projections and histories

###  create lists of historical runs needed for each RCP
lf26 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/26/')
lf26 <- gsub('rcp26', 'historical', lf26)
lf26 <- gsub('200601-210012', '*', lf26)
lf26 <- gsub('200601-210011', '*', lf26)
lf26 <- gsub('200601-209912', '*', lf26)
lf26 <- gsub('200512-209912', '*', lf26)
lf26 <- gsub('200512-209911', '*', lf26)
lf26 <- gsub('200512-210011', '*', lf26)
lf26 <- gsub('200512-210012', '*', lf26)

lf45 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/45/')
lf45 <- gsub('rcp45', 'historical', lf45)
lf45 <- gsub('200601-210012', '*', lf45)
lf45 <- gsub('200601-210011', '*', lf45)
lf45 <- gsub('200601-209912', '*', lf45)
lf45 <- gsub('200512-209912', '*', lf45)
lf45 <- gsub('200512-209911', '*', lf45)
lf45 <- gsub('200512-210011', '*', lf45)
lf45 <- gsub('200512-210012', '*', lf45)

lf60 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/60/')
lf60 <- gsub('rcp60', 'historical', lf60)
lf60 <- gsub('200601-210012', '*', lf60)
lf60 <- gsub('200601-210011', '*', lf60)
lf60 <- gsub('200601-209912', '*', lf60)
lf60 <- gsub('200512-209912', '*', lf60)
lf60 <- gsub('200512-209911', '*', lf60)
lf60 <- gsub('200512-210011', '*', lf60)
lf60 <- gsub('200512-210012', '*', lf60)

lf85 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/85/')
lf85 <- gsub('rcp85', 'historical', lf85)
lf85 <- gsub('200601-210012', '*', lf85)
lf85 <- gsub('200601-210011', '*', lf85)
lf85 <- gsub('200601-209912', '*', lf85)
lf85 <- gsub('200512-209912', '*', lf85)
lf85 <- gsub('200512-209911', '*', lf85)
lf85 <- gsub('200512-210011', '*', lf85)
lf85 <- gsub('200512-210012', '*', lf85)

lf26.DJF <- lf26[grep('DJF', lf26)]
lf26.JJA <- lf26[grep('JJA', lf26)]
lf45.DJF <- lf45[grep('DJF', lf45)]
lf45.JJA <- lf45[grep('JJA', lf45)]
lf60.DJF <- lf60[grep('DJF', lf60)]
lf60.JJA <- lf60[grep('JJA', lf60)]
lf85.DJF <- lf85[grep('DJF', lf85)]
lf85.JJA <- lf85[grep('JJA', lf85)]


write.csv(lf26.DJF, row.names=F)
write.csv(lf45.DJF, row.names=F)
write.csv(lf60.DJF, row.names=F)
write.csv(lf85.DJF, row.names=F)
write.csv(lf26.JJA, row.names=F)
write.csv(lf45.JJA, row.names=F)
write.csv(lf60.JJA, row.names=F)
write.csv(lf85.JJA, row.names=F)

#### create lists of files for DJF and JJA
lf26 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/26/')
lf26.DJF <- lf26[grep('DJF', lf26)]
lf26.JJA <- lf26[grep('JJA', lf26)]

lf45 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/45/')
lf45.DJF <- lf45[grep('DJF', lf45)]
lf45.JJA <- lf45[grep('JJA', lf45)]

lf60 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/60/')
lf60.DJF <- lf60[grep('DJF', lf60)]
lf60.JJA <- lf60[grep('JJA', lf60)]

lf85 <- list.files('/Volumes/CCDatStat/CMIP5/global_mon/BCSD/upshot/85/')
lf85.DJF <- lf85[grep('DJF', lf85)]
lf85.JJA <- lf85[grep('JJA', lf85)]

write.csv(lf26.DJF, row.names=F)
write.csv(lf45.DJF, row.names=F)
write.csv(lf60.DJF, row.names=F)
write.csv(lf85.DJF, row.names=F)
write.csv(lf26.JJA, row.names=F)
write.csv(lf45.JJA, row.names=F)
write.csv(lf60.JJA, row.names=F)
write.csv(lf85.JJA, row.names=F)
