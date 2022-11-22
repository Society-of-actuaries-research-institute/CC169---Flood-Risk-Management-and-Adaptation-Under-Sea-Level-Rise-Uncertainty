
load.lib<-c("zoo","lubridate","zoo","xtable","dotCall64","rlang","dplyr","moments","zeallot",
            "TideHarmonics", "forecast", 'tseries','extRemes','quantmod')
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)


