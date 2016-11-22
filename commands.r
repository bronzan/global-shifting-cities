lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
lonlat <- expand.grid(lon, lat)
tmax.df <- data.frame(lonlat)
tmax.array <- ncvar_get(nc, "tasmax")
fillvalue <- 1.00000002004088e+20;
tmax.array[tmax.array == fillvalue] <- NA;
slice <- as.vector(tmax.array[,,433])
tmax.df <- cbind(tmax.df, slice)
startpt <- gregexpr("_", fn)[[1]][4]
endpt <- gregexpr("_", fn)[[1]][5]
substr(fn, startpt + 1, endpt - 1)
