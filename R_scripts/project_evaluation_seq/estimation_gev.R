# Obtain annual maximum
data <- read.csv(file = paste(opath,"hourly_series_tide", ".csv", sep=""), sep=",", header=TRUE) %>%
  mutate(time=as.POSIXct(time,format = "%Y-%m-%d %H:%M",tz="UTC")) %>%
  mutate(year = year(time)) %>%
  group_by(year) %>% 
  summarise(wldsm = max(wlds), msl = mean(wlds))%>%
  mutate(wlds = wldsm - msl)%>%
  filter(year> 1979)

# add the most recent mean sea level back to deseaonalized, demeaned water level
data$wlds = data$wlds + data$msl[nrow(data)]
# ADF test
plot(data$year, data$wlds)
lines(data$year, data$wlds)
plot(data$year, data$msl)
adf.test(data$wlds)

# GEV estimation
#=== Plot sea level ===
png(file = paste(fig, "qldsl.png", sep=""), res = 200, width = 1200, height = 800)
plot(c(min(data$year), max(data$year)), c(min(data[,2:3]), max(data[,2:3])),xlab = 'Year', ylab = 'Sea level (mm)', cex=0)
points(data$year, data$msl, col='black')
lines(data$year, data$msl, col='black')
points(data$year, data$wldsm, col='blue', lty=1.5)
lines(data$year, data$wldsm, col='blue', lty=1.5)  # wldsm is annual maximum that is deseasonalized.
dev.off()

#=== Estimate GEV distribution ===
# GEV for excess of max above mean sea level 
fit <- fevd(data$wlds, type = 'GEV', units = "MM")
fit
plot(fit)
cie0 <- ci(fit, return.period = c(10,100, 1000, 10000))
cie0
#=== Brownian motion for mean sea level ===
y <- data$msl[2:nrow(data)] - data$msl[1:(nrow(data)-1)]
mean(y)
sd(y)
