timesq <- function(date, halfday){
  #timesq provides hourly sequence of time in a day and in a night.
  if(halfday=='1'){
    starthour <- as.POSIXct(paste(date,'00:00'),format = "%Y-%m-%d %H:%M",  tz="UTC")
    endhour <- as.POSIXct(paste(date,'11:00'), format = "%Y-%m-%d %H:%M", tz="UTC")
    out <- seq(starthour, endhour, by='hours')
  }
  if(halfday=='2'){
    starthour <- as.POSIXct(paste(date,'12:00'), format = "%Y-%m-%d %H:%M", tz="UTC")
    endhour <- as.POSIXct(paste(date,'23:00'), format = "%Y-%m-%d %H:%M",  tz="UTC")
    out <- seq(starthour, endhour, by='hours')
  }
  return(out)
}