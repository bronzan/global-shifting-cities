state.trends <- read.csv("~/Desktop/Tables/seasonal trends/state_trends.csv")
state.meta <- read.csv("~/Desktop/Tables/seasonal trends/state_meta.csv")
state.data <- read.csv("~/Desktop/Tables/seasonal trends/state_data.csv")


# class(state.trends$State) <- 'integer'
# class(state.meta$State) <- 'integer'

state.trends <- merge(state.trends, state.meta)


winter.blue <- "#00c8ff"
spring.green <- "#8fdf26"
summer.red <- "#e02d1f"
fall.orange <- "#fba924"

xvals <- 0:44
szns <- unique(state.trends$Season)
stateNums <- unique(state.trends$State)
slp <- c()
yvals <- c()
#pdf("~/Desktop/Plots/Seasonal State Warming Trends 2015-01-22.pdf", width=8, height=6)
	for ( st in stateNums ) {
		stateName <- unique(subset(state.trends, State==st)$StateName)
		#pdf(paste0("~/Desktop/Plots/Seasonal State Warming Trends/Seasonal_State_Trends_PMA_out/Seasons-2015-12-11/",stateName,".pdf"), width=10, height=6)
		pdf(paste0("~/Desktop/Plots/Seasonal State Warming Trends/Seasonal_State_Trends_PMA_out/Grids-2015-12-11/",stateName,".pdf"), width=10, height=6)
			rm(yvals)
			yvals <- c()
			for ( szn in szns ) {
				rm(slp)
				slp <- subset(state.trends, State==st & Season==szn)$Slope[1]
				yvals <- cbind(yvals, slp*xvals)
			}
			min.slp <- min(c(0, yvals[45,]))
			max.slp <- max(yvals[45,])
			colnames(yvals) <- szns
	#		plot(xvals, yvals[,1], type = 'n', ylim=c(min(yvals),max(yvals)), main = paste(stateName, 'trends'))
#			plot(xvals, yvals[,1], xlab="Number of Years since 1970", ylab="Average Temperature Increase since 1970", type = 'n', ylim=c(-0.5,4.0), main = paste(stateName, 'trends'))
#			plot(xvals, yvals[,1], xlab='', ylab='', type = 'n', ylim=c(-0.5,4.0), axes=F, bty='n')
			plot(xvals, yvals[,1], xlab='', ylab='', type = 'n', ylim=c(min.slp, max.slp), axes=F, bty='n')
# 			axis(2, las=1, col.axis="#cacaca", font=2, tcl=0, lty='blank')
# 			grid(nx=NA, ny=NULL, lwd=4, lty="11")
# 			winter blue
# 			spring green
# 			summer red
# 			fall orange
			lines(xvals, yvals[,"SON"], col=fall.orange, lwd=8)
			lines(xvals, yvals[,"JJA"], col=summer.red, lwd=8)
			lines(xvals, yvals[,"MAM"], col=spring.green, lwd=8)
			lines(xvals, yvals[,"DJF"], col=winter.blue, lwd=8)
#			legend(0, 3.5, c("Winter", "Spring", "Summer", "Fall"), lwd=4, col = c(winter.blue, spring.green, summer.red, fall.orange))
		dev.off()
	}
#dev.off()


### For interactive...

steepest <- reshape(state.trends, v.names = c("Slope","Intercept","P_value"), idvar = "StateAbb", timevar = "Season", direction = "wide")
steepest$littleAbb <- tolower(steepest$StateAbb)
steepest$steepest <- NA
for (i in 1:nrow(steepest)) {
	steepest$steepest[i] <- which.max(c(steepest$Slope.DJF[i], steepest$Slope.MAM[i],
										steepest$Slope.JJA[i], steepest$Slope.SON[i]))
}
write.csv(steepest, "~/Desktop/Tables/steepest_seasonal_trends.csv")

# Shared https://docs.google.com/a/climatecentral.org/spreadsheets/d/1DENPCNumWl3m9JQqrAzeA8GMgGOKYC7BByLfPyuzKUE/edit?usp=sharing
