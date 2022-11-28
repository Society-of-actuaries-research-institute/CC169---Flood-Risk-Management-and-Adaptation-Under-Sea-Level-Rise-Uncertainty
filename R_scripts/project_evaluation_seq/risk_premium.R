sm <- read.table(paste(opath, 'allord.csv', sep=""), sep=",", header=TRUE)%>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"))%>%
  mutate(year = year(date))%>%
  group_by(year)%>%
  slice(which.max(date))%>%
  ungroup()%>%
  select(date, price)%>%
  mutate(smr = (price - Lag(price,1))/Lag(price,1),year=year(date))

sl <- read.csv(file = paste(opath,"hourly_series_tide", ".csv", sep=""), sep=",", header=TRUE)%>%
  mutate(time=as.POSIXct(time,format = "%Y-%m-%d %H:%M",tz="UTC")) %>%
  mutate(year = year(time)) %>%
  filter(year < 2019, year>1945 ) %>%
  group_by(year) %>% 
  summarise(msl = mean(wlds))%>%
  select(year, msl)%>%
  mutate(slr = msl - Lag(msl,1))%>%
  inner_join(sm,by='year')%>%
  na.omit()
cor(sl$slr, sl$smr)
cor(sl$slr, sl$smr)*0.1274*31.29  #market risk premium