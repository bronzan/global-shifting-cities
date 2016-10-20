## Unpacks nc files in "data" to extract just 1986-2005 (monthly data)

library(ncdf4)
print(Sys.time())



options(stringsAsFactors=FALSE)
fns <- list.files("data")


## generate column names
yrs <- 1986:2005
mos <- c("01","02","03","04","05","06","07","08","09","10","11","12")
col_labels <- c()
for (y in 1:length(yrs)) {
	col_labels <- c(col_labels, paste0(yrs[y], mos))
}



## Loop through all files in directory, unpacking each
for (f in 1:length(fns)) {

	##filename to open
	fn <- paste0("data/", fns[f])

	##get model name
	startpt <- gregexpr("_", fn)[[1]][4]
	endpt <- gregexpr("_", fn)[[1]][5]
	model_name <- substr(fn, startpt + 1, endpt - 1)

	## Open netcdf, get lat and lon and prepare data frame
	nc <- nc_open(fn)
	lon <- ncvar_get(nc, "longitude")
	lat <- ncvar_get(nc, "latitude")
	lonlat <- expand.grid(lon, lat)
	tmax.df <- data.frame(lonlat)

	## Get 
	tmax.array <- ncvar_get(nc, "tasmax")
	fillvalue <- 1.00000002004088e+20;
	tmax.array[tmax.array == fillvalue] <- NA;
	for (s in 433:672) {
		slice <- as.vector(tmax.array[,,433])
		tmax.df <- cbind(tmax.df, slice)
		colnames(tmax.df)[ncol(tmax.df)] <- col_labels[s - 432]
	}

	save(tmax.df, file="tmx_monthly_1986_to_2005_", modelname, ".RData")
}



print(Sys.time())