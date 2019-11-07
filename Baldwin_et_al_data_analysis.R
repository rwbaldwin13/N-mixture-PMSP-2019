# Comparison of Analysis Methods for White-tailed Deer Population Surveys
# Baldwin et al. 2019
# Prepared November 2019


# ----------------------------------------- FLIR VALIDATION DENSITY ESTIMATION --------------------------------------------------- #

# remove old objects
  rm(list=ls())

# load library
  library(car)

# read-in observer data
  flir.dat <- read.csv("pub_data/flir_obs.csv")
  transect.dat <- read.csv("pub_data/transect_information.csv")
  
# set seed for reproducibility (applies for all future analyses as well)
  set.seed(4)

# merge dataframes for anaysis
  flir.dat <- merge(flir.dat, transect.dat, by.x = "TRANSECT", by.y = "TRANSECT")

# convert counts to densities
  flir.dat$DENSITY <- (flir.dat$COUNT / flir.dat$AREA)

# subset data by flight
  f.1 <- subset(flir.dat, FLIGHT == 1, drop = T)
  f.2 <- subset(flir.dat, FLIGHT == 2, drop = T)
  f.3 <- subset(flir.dat, FLIGHT == 3, drop = T)
  f.4 <- subset(flir.dat, FLIGHT == 4, drop = T)
  f.5 <- subset(flir.dat, FLIGHT == 5, drop = T)

# aggregate data by transect and observer
  f1trans <- aggregate(f.1$DENSITY, by = list(f.1$TRANSECT, f.1$OBSERVER), FUN = "sum")
  names(f1trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f1trans <- aggregate(f1trans$DENSITY, by = list(f1trans$TRANSECT), FUN = 'mean')
  names(f1trans) <- c("TRANSECT","DENSITY")

  f2trans <- aggregate(f.2$DENSITY, by = list(f.2$TRANSECT, f.2$OBSERVER), FUN = "sum")
  names(f2trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f2trans <- aggregate(f2trans$DENSITY, by = list(f2trans$TRANSECT), FUN = 'mean')
  names(f2trans) <- c("TRANSECT","DENSITY")
  
  f3trans <- aggregate(f.3$DENSITY, by = list(f.3$TRANSECT, f.3$OBSERVER), FUN = "sum")
  names(f3trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f3trans <- aggregate(f3trans$DENSITY, by = list(f3trans$TRANSECT), FUN = 'mean')
  names(f3trans) <- c("TRANSECT","DENSITY")
  
  f4trans <- aggregate(f.4$DENSITY, by = list(f.4$TRANSECT, f.4$OBSERVER), FUN = "sum")
  names(f4trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f4trans <- aggregate(f4trans$DENSITY, by = list(f4trans$TRANSECT), FUN = 'mean')
  names(f4trans) <- c("TRANSECT","DENSITY")

  f5trans <- aggregate(f.5$DENSITY, by = list(f.5$TRANSECT, f.5$OBSERVER), FUN = "sum")
  names(f5trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f5trans <- aggregate(f5trans$DENSITY, by = list(f5trans$TRANSECT), FUN = 'mean')
  names(f5trans) <- c("TRANSECT","DENSITY")

# join data to calculate average density between flights
  density.dat <- rbind(f1trans, f2trans, f3trans, f4trans, f5trans)
  mean(density.dat$DENSITY)
  density <- c(mean(f1trans$DENSITY), mean(f2trans$DENSITY), mean(f3trans$DENSITY), mean(f4trans$DENSITY), mean(f5trans$DENSITY))

# create dataframe for modelling
  fl.agg <- data.frame(Flight = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6)), 
                       Transect = c(f1trans$TRANSECT, f2trans$TRANSECT, f3trans$TRANSECT, f4trans$DENSITY, f5trans$DENSITY), 
                       Density = c(f1trans$DENSITY, f2trans$DENSITY, f3trans$DENSITY, f4trans$DENSITY, f5trans$DENSITY))
  
# bootstrapped CI for density
  bstrap <- NULL
  for (i in 1:10000){
  bstrap <- c(bstrap, mean(sample(density, 5,replace=T)))}
  rbind(c(quantile(bstrap,.025), quantile(bstrap, .975)))


# run ANOVA on flight estimates
  m1 <- lm(fl.agg$Density ~ fl.agg$Flight)
  Anova(m1, typ = '2')

# ----------------------------------------- MR BOOTSTRAPPING FOR PRECISON ESTIMATION --------------------------------------------------- #

# remove old objects
  rm(list = ls())
  
# read in data (UM stands for number of unique males at each site)
  mr.dat <- read.csv("pub_data/mr_error.csv")
  head(mr.dat)
  str(mr.dat)
  
# simulate survey 10,000 times 
  set.seed(4)
  
  N <- NULL
  sp <- NULL
  for(i in 1:1000){
    sp <- data.frame(SITE = sample(c(1:22), 22, replace = T))
    mg <- merge(sp, mr.dat, by.y = 'SITE')
    correction <- sum(mg$MALE) / sum(mg$UM)
    FC <- sum(mg$FEMALE)/correction
    YC <- sum(mg$FAWN)/correction
    N[i] <- sum(c(FC, YC, sum(mg$UM)))
  }
  
  quantile(N, 0.025) / 10.24
  quantile(N, 0.975) / 10.24
  sd(N)
  


#  ----------------------------------------- N-MIXTURE MODELLING DEER ABUNDANCE JAN - MAR 2018 --------------------------------------------------- #
    
# remove old objects
  rm(list=ls())
  
# load libraries
  library(unmarked)
  library(AICcmodavg)
  library(reshape)

# read-in data
  deer.dat <- read.csv("pub_data/deer_obs_jan_march.csv")
  cam.dat <- read.csv("pub_data/camera_covariates.csv")

# create list of ordered levels for sorting
  out <- NULL
  for(i in 1:22) {
  tmp <- noquote(paste("PM",i,sep=""))
  out <- append(out, tmp)
  } 
  deer.dat$SITE <- ordered(deer.dat$SITE, levels=out)

# convert DATE.TIME to POSIX
  deer.dat$DATE.TIME <- as.POSIXct(strptime(deer.dat$DATE.TIME, format = "%m/%d/%Y %H:%M"))

# order observations by site and time
  deer.dat <- deer.dat[order(deer.dat[,2], deer.dat[,9]),]

# aggregate by max count by hour prior to merge
  head(deer.dat, 50)
  deer.dat$DATE.TIME.HOUR <- round.POSIXt(deer.dat$DATE.TIME, "hours")
  count.by.hour <- aggregate(deer.dat$COUNT, by = list(deer.dat$SITE, as.character(deer.dat$DATE.TIME.HOUR)), FUN = max)
  names(count.by.hour) <- c("SITE", "DATE.HOUR", "COUNT")
  head(count.by.hour)

# make continuous calender of all observation periods for merge
  all.sites <- unique(deer.dat$SITE)

# the survey began on Jan 1st at midnight, thus the start date
  min.date.time <- as.POSIXct(strptime("2018-01-01 00:00:00", "%Y-%m-%d %H:%M:%S"))
  max.date.time <- as.POSIXct(strptime("2018-04-01 00:00:00", "%Y-%m-%d %H:%M:%S"))           
  seq.date.time <- seq.POSIXt(from = min.date.time, to = max.date.time, by = "hour")
  rep.date.time <- rep(seq.date.time, length(all.sites))

# create site reps for new data frame
  rep.sites <- sort(rep(all.sites, length(seq.date.time)))

# create data frame for merge
  survey.sites.df <- data.frame(SITE = rep.sites, DATE.HOUR = rep.date.time)
  head(survey.sites.df)

# merge site by survey df with aggregated count data
  deer.df <- merge(x = survey.sites.df, y = count.by.hour, by.x = c("SITE","DATE.HOUR"), by.y = c("SITE","DATE.HOUR"), all.x = T, all.y = T)
  deer.df[is.na(deer.df)] <- 0
  
# trim dataset so it is exactly 89 days
  days.on <- as.numeric(round(max.date.time - min.date.time, 0))
  temp <- deer.df[order(deer.df$DATE.HOUR),]
  rownames(temp) <- 1:nrow(temp)
  temp <- temp[1:47520,]
  temp <- temp[order(temp$SITE, temp$DATE.HOUR),]
  rownames(temp) <- 1:nrow(temp)
  temp$DATE <- as.Date.POSIXct(temp$DATE.HOUR, format = "%Y-%m-%d")

# make duplicate dataframes for sensitivity analysis of survey length
  temp.1 <- temp
  temp.2 <- temp
  temp.3 <- temp
  temp.4 <- temp
  temp.6 <- temp
  temp.8 <- temp
  temp.12 <- temp
  temp.day <- temp
  temp.48 <- temp
  temp.5day <- temp
  temp.week <- temp


# add new column to aggregate data for modelling
  temp.week$AGG <-sort(rep(1:12, 22*7.5))     # 12, 7-day survey blocks
  temp.5day$AGG <- sort(rep(1:18, 120))       # 18, 5-day survey blocks
  temp.48$AGG <- sort(rep(1:45, 48))          # 45, 2-day survey blocks
  temp.day$AGG <- sort(rep(1:days.on, 24))    # 90, 1 day survey blocks
  temp.12$AGG <- rep(sort(rep(1:2, 12)), 22*90)# 2, 12-hour survey blocks per day
  temp.8$AGG <- rep(sort(rep(1:3, 8)), 22*90) # 3, 8-hour survey blocks per day
  temp.6$AGG <- rep(sort(rep(1:4, 6)), 22*90) # 4, 6-hour survey blocks per day
  temp.4$AGG <- rep(sort(rep(1:6, 4)), 22*90) # 6, 4-hour survey blocks per day
  temp.3$AGG <- rep(sort(rep(1:8, 3)), 22*90) # 8, 3-hour survey blocks per day
  temp.2$AGG <- rep(sort(rep(1:12, 2)), 22*90) # 12, 2-hour blocks per day
  temp.1$AGG <- rep(sort(rep(1:24, 1)), length(all.sites)*days.on) #24, 1-hour blocks per day
  temp$AGG <- rep(sort(rep(1:24, 1)), length(all.sites)*days.on) # 24, 1 hour blocks per day for model selection and comparisons

# aggregate data based #on AGG column
  temp.agg <- aggregate(COUNT ~ SITE  + AGG + DATE, data = temp, FUN = "max")
  temp.1agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.1, FUN = "max")
  temp.2agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.2, FUN = "max")
  temp.3agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.3, FUN = "max")
  temp.4agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.4, FUN = "max")
  temp.6agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.6, FUN = "max")
  temp.8agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.8, FUN = "max")
  temp.12agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.12, FUN = "max")
  temp.dayagg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.day, FUN = "max")
  temp.48agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.48, FUN = "max")
  temp.5dayagg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.5day, FUN = "max")
  temp.weekagg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.week, FUN = "max")
  deer.df <- temp.agg

# create new column which is a combination of DATE and AGG
  temp.1agg$DATE.AGG <- paste(temp.1agg$DATE, temp.1agg$AGG, sep = '.')
  temp.2agg$DATE.AGG <- paste(temp.2agg$DATE, temp.2agg$AGG, sep = '.')
  temp.3agg$DATE.AGG <- paste(temp.3agg$DATE, temp.3agg$AGG, sep = '.')
  temp.4agg$DATE.AGG <- paste(temp.4agg$DATE, temp.4agg$AGG, sep = '.')
  temp.6agg$DATE.AGG <- paste(temp.6agg$DATE, temp.6agg$AGG, sep = '.')
  temp.8agg$DATE.AGG <- paste(temp.8agg$DATE, temp.8agg$AGG, sep = '.')
  temp.12agg$DATE.AGG <- paste(temp.12agg$DATE, temp.12agg$AGG, sep = '.')
  temp.dayagg$DATE.AGG <- paste(temp.dayagg$DATE, temp.dayagg$AGG, sep = '.')
  temp.48agg$DATE.AGG <- paste(temp.48agg$DATE, temp.48agg$AGG, sep = '.')
  temp.5dayagg$DATE.AGG <- paste(temp.5dayagg$DATE, temp.5dayagg$AGG, sep = '.')
  temp.weekagg$DATE.AGG <- paste(temp.weekagg$DATE, temp.weekagg$AGG, sep = '.')
  deer.df$DATE.AGG <- paste(deer.df$DATE, deer.df$AGG, sep=".")

# create simplified df of counts
  count.1 <- temp.1agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.2 <- temp.2agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.3 <- temp.3agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.4 <- temp.4agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.6 <- temp.6agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.8 <- temp.8agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.12 <- temp.12agg[, c("SITE", 'DATE.AGG', 'COUNT')]
  count.day <- temp.dayagg[,c("SITE", 'DATE.AGG', 'COUNT')]
  count.48 <- temp.48agg[,c("SITE", 'DATE.AGG', 'COUNT')]
  count.5day <- temp.5dayagg[,c("SITE", 'DATE.AGG', 'COUNT')]
  count.week <- temp.weekagg[,c("SITE", 'DATE.AGG', 'COUNT')]

# cast data into matrices for unmarked package
  count.1 <- cast(data = count.1, formula = SITE~DATE.AGG,mean)
  count.2 <- cast(data = count.2, formula = SITE~DATE.AGG, mean)
  count.3 <- cast(data = count.3, formula = SITE~DATE.AGG, mean)
  count.4 <- cast(data = count.4, formula = SITE~DATE.AGG, mean)
  count.6 <- cast(data = count.6, formula = SITE~DATE.AGG, mean)
  count.8 <- cast(data = count.8, formula = SITE~DATE.AGG, mean)
  count.12 <- cast(data = count.12, formula = SITE~DATE.AGG, mean)
  count.day <- cast(data = count.day, formula = SITE~DATE.AGG)
  count.48 <- cast(data = count.48, formula = SITE~DATE.AGG)
  count.5day <- cast(data = count.5day, formula = SITE~DATE.AGG)
  count.week <- cast(data = count.week, formula = SITE~DATE.AGG, mean)

  count.1 <- as.matrix.cast_df(count.1,-1)
  count.2 <- as.matrix.cast_df(count.2,-1)
  count.3 <- as.matrix.cast_df(count.3,-1)
  count.4 <- as.matrix.cast_df(count.4,-1)
  count.6 <- as.matrix.cast_df(count.6,-1)
  count.8 <- as.matrix.cast_df(count.8,-1)
  count.12 <- as.matrix.cast_df(count.12,-1)
  count.day <- as.matrix.cast_df(count.day,-1)
  count.48 <- as.matrix.cast_df(count.48,-1)
  count.5day <- as.matrix.cast_df(count.5day,-1)
  count.week <- as.matrix.cast_df(count.week,-1)

# repeat process for model selection
  count.df <- deer.df[, c("SITE","DATE.AGG","COUNT")]
  count.df <- cast(data = count.df, formula = SITE~DATE.AGG, mean) #timemeans <- cast(mdata, time~variable, mean) 
  count.df <- as.matrix.cast_df(count.df[,-1]) 

# create site covariate dataframe
  deer.covars <- data.frame(ELE=scale(cam.dat$site.ele), ASPECT = cam.dat$slopecor, SLOPE = scale(cam.dat$slope), EDGE=scale(cam.dat$EDGE), ROAD=scale(cam.dat$ROAD))
  str(deer.covars)
  head(deer.covars)

# create unmarked frame with data 
  pcount.list <- list(count.1, count.2, count.3, count.4, count.6, count.8, count.12, count.day, count.48, count.5day, count.week)
  pcount.list <- lapply(pcount.list, function(x) unmarkedFramePCount(y = x, siteCovs = deer.covars))
  deer.pcount <- unmarkedFramePCount(y = count.df, siteCovs = deer.covars)
  deer.day <- unmarkedFramePCount(y = count.day, siteCovs = deer.covars)

# model of Royal 2004b based on counts
# general formula is: model <- pcount(~detection_formula ~occupancy_formula, dataframe, K=100, se=TRUE)
# starting values generated by iteratively running models

  m0 <-  pcount(~1 ~1, data = deer.pcount, K = 150, se = T, mixture = "P",
              starts = c(2.53, -6.26),
              control = list(trace = T, REPORT = 1))


  m1 <-  pcount(~1 ~ROAD, data = deer.pcount, K = 150, se = T, 
              starts = c(2.53, 0,-6.26),
              control = list(trace=T, REPORT = 1)) 

  m2 <- pcount(~1 ~EDGE, data = deer.pcount, K = 150, se = T, 
             starts = c(2.53, 0,-6.26), 
             control = list(trace=T, REPORT = 1))

  m3 <- pcount(~1 ~ ELE, data = deer.pcount, K = 150, se = T, 
             starts = c(2.53, 0, -6.26),
             control = list(trace=T, REPORT = 1)) 

  m4 <- pcount(~1 ~ASPECT, data = deer.pcount, K = 150, se = T, 
             starts = c(2.4, 0.3, -6.26), 
             control = list(trace=T, REPORT = 1)) 

  m5 <- pcount(~1 ~ASPECT + EDGE, data = deer.pcount, K = 150, se = T, 
             starts = c(2.4, 0.3, 0,-6.26), 
             control = list(trace=T, REPORT = 1))

  m6 <- pcount(~1 ~ASPECT*EDGE, data = deer.pcount, K = 150, se = T, 
             starts = c(2.4,  0.3, 0, 0,-6.26), 
             control = list(trace=T, REPORT = 1))

  m7 <- pcount(~1 ~ASPECT*ELE, data = deer.pcount, K = 150, se = T,
             starts = c(2.4, 0.3, 0, 0, -6.26),
             control = list(trace=T, REPORT = 1)) 

  m8 <- pcount(~SLOPE ~1, deer.pcount, K = 150, se = T, 
             starts = c(2.53, -6.26, 0), 
             control = list(trace=T, REPORT = 1)) 

  m9<- pcount(~ASPECT*SLOPE ~1, deer.pcount, K = 150, se = T, 
            starts = c(2.4,-6.26, 0.26, 0, 0.05), 
            control = list(trace=T, REPORT = 1))


  m10 <- pcount(~ASPECT*SLOPE ~ASPECT*EDGE, deer.pcount, K = 150, se = T, 
              starts = c(2.4, 0.26, 0, 0.05, -6.26, -0.7,0, 0), 
              control = list(trace=T, REPORT = 1))

  m11 <- pcount(~ASPECT*SLOPE ~EDGE, data = deer.pcount, K = 150, se = T, 
              starts = c(2.5, 0, -5.7, -0.7, 0, 0), 
              control = list(trace=T, REPORT = 1))


  m12<- pcount(~ASPECT*SLOPE ~ELE, deer.pcount, K = 150, se = T, 
             starts = c(2.53, 0, -5.72, -0.7, 0, 0 ), 
             control = list(trace=T, REPORT = 1))

  m13<- pcount(~ASPECT*SLOPE ~ASPECT, deer.pcount, K = 150, se = T, 
             starts = c(2.4,0.26, -5.72, -0.72, -0.05, 0.07), 
             control = list(trace=T, REPORT = 1))


# create fit list of models
  fl <- modSel(fitList(Null = m0, A = m1,B = m2, C = m3, D = m4, E = m5, G = m6, H = m7, I = m8,K= m9,L = m10, M = m11, N = m12, O = m13))
  
# calculate aics of model and select top model
  fl
  summary(m13)

# run null model on different aggregations as sensitivity analysis for survey length
  model.list <- lapply(pcount.list, function(x) pcount(~1 ~1, x, K = 150, se = T, 
                                                     starts = c(2.5,-6.26), 
                                                     control = list(trace=T, REPORT = 1)) )  


# calculate final estimates for comparison between methods
  re <- ranef(m13)
  EBUP <- bup(re, stat="mean")
  CI<- confint(re, level=0.95)
  Nmixture.estimate <- rbind(c(Estimate = sum(EBUP), colSums(CI)))
  Nmixture.estimate /10.24


# repeat process on the sensitivity analysis models
  re.list <- lapply(model.list, function(x) ranef(x))
  EBUP.list <- lapply(re.list, function (x) bup(x, stat = 'mean'))
  CI.list <- lapply(re.list, function(x) confint(x, level = 0.95))

  hour1.EBUP <- unlist(EBUP.list[[1]])
  hour1.CI <- unlist(CI.list[[1]])
  hour1.CI <- as.data.frame(hour1.CI)
  hour1.EST <- data.frame(Estimate = sum(hour1.EBUP), Lower = sum(hour1.CI$`2.5%`), Upper = sum(hour1.CI$`97.5%`))


  hour2.EBUP <- unlist(EBUP.list[[2]])
  hour2.CI <- unlist(CI.list[[2]])
  hour2.CI <- as.data.frame(hour2.CI)
  hour2.EST <- data.frame(Estimate = sum(hour2.EBUP), Lower = sum(hour2.CI$`2.5%`), Upper = sum(hour2.CI$`97.5%`))

  hour3.EBUP <- unlist(EBUP.list[[3]])
  hour3.CI <- unlist(CI.list[[3]])
  hour3.CI <- as.data.frame(hour3.CI)
  hour3.EST <- data.frame(Estimate = sum(hour3.EBUP), Lower = sum(hour3.CI$`2.5%`), Upper = sum(hour3.CI$`97.5%`))

  hour4.EBUP <- unlist(EBUP.list[[4]])
  hour4.CI <- unlist(CI.list[[4]])
  hour4.CI <- as.data.frame(hour4.CI)
  hour4.EST <- data.frame(Estimate = sum(hour4.EBUP), Lower = sum(hour4.CI$`2.5%`), Upper = sum(hour4.CI$`97.5%`))

  hour6.EBUP <- unlist(EBUP.list[[5]])
  hour6.CI <- unlist(CI.list[[5]])
  hour6.CI <- as.data.frame(hour6.CI)
  hour6.EST <- data.frame(Estimate = sum(hour6.EBUP), Lower = sum(hour6.CI$`2.5%`), Upper = sum(hour6.CI$`97.5%`))

  hour8.EBUP <- unlist(EBUP.list[[6]])
  hour8.CI <- unlist(CI.list[[6]])
  hour8.CI <- as.data.frame(hour8.CI)
  hour8.EST <- data.frame(Estimate = sum(hour8.EBUP), Lower = sum(hour8.CI$`2.5%`), Upper = sum(hour8.CI$`97.5%`))

  hour12.EBUP <- unlist(EBUP.list[[7]])
  hour12.CI <- unlist(CI.list[[7]])
  hour12.CI <- as.data.frame(hour12.CI)
  hour12.EST <- data.frame(Estimate = sum(hour12.EBUP), Lower = sum(hour12.CI$`2.5%`), Upper = sum(hour12.CI$`97.5%`))

  day.EBUP <- unlist(EBUP.list[[8]])
  day.CI <- unlist(CI.list[[8]])
  day.CI <- as.data.frame(day.CI)
  day.EST <- data.frame(Estimate = sum(day.EBUP), Lower = sum(day.CI$`2.5%`), Upper = sum(day.CI$`97.5%`))

  hour48.EBUP <- unlist(EBUP.list[[9]])
  hour48.CI <- unlist(CI.list[[9]])
  hour48.CI <- as.data.frame(day.CI)
  hour48.EST <- data.frame(Estimate = sum(hour48.EBUP), Lower = sum(hour48.CI$`2.5%`), Upper = sum(hour48.CI$`97.5%`))

  day5.EBUP <- unlist(EBUP.list[[10]])
  day5.CI <- unlist(CI.list[[10]])
  day5.CI <- as.data.frame(day.CI)
  day5.EST <- data.frame(Estimate = sum(day5.EBUP), Lower = sum(day5.CI$`2.5%`), Upper = sum(day5.CI$`97.5%`))

  week.EBUP <- unlist(EBUP.list[[11]])
  week.CI <- unlist(CI.list[[11]])
  week.CI <- as.data.frame(day.CI)
  week.EST <- data.frame(Estimate = sum(week.EBUP), Lower = sum(week.CI$`2.5%`), Upper = sum(week.CI$`97.5%`))


  window.est <- rbind(hour1.EST, hour2.EST, hour3.EST, hour4.EST, hour6.EST, hour8.EST, hour12.EST, day.EST, hour48.EST, day5.EST, week.EST)
  window.est$Window <- c('1 hour','2 hour','3 hour', '4 hour', '6 hour', '8 hour', '12 hour', '24 hour', '48hour', '5day', '1 week')


#  ----------------------------------------- N-MIXTURE MODELLING DEER ABUNDANCE SUMMER16 - SPRING 2018 --------------------------------------------------- #
# remove old objects
  rm(list=ls())
  
# read in new data
  deer.dat <- read.csv("deer_obs_full_period.csv")
  cam.dat <- read.csv("camera_covariates.csv")
  
# create list of ordered levels for sorting later
  out <- NULL
  for(i in 1:22) {
    tmp <- noquote(paste("PM",i,sep=""))
    out <- append(out, tmp)
  } 
  deer.dat$SITE <- ordered(deer.dat$SITE, levels=out)
  
# convert date to POSIXct
  deer.dat$DATE.TIME <- as.POSIXct(deer.dat$DATE.TIME, format = '%m/%d/%Y %H:%M:%S')
  
# order observations by site and time
  deer.dat <- deer.dat[order(deer.dat$SITE, deer.dat$DATE.TIME),]

# aggregate by max count by hour prior to merge
  head(deer.dat, 50)
  deer.dat$DATE.TIME.HOUR <- round.POSIXt(deer.dat$DATE.TIME, "hours")
  count.by.hour <- aggregate(deer.dat$COUNT, by = list(deer.dat$SITE, as.character(deer.dat$DATE.TIME.HOUR)), FUN = max)
  names(count.by.hour) <- c("SITE", "DATE.HOUR", "COUNT")
  head(count.by.hour)
  count.by.hour <- count.by.hour[order(count.by.hour$SITE, count.by.hour$DATE.HOUR),]
  
# make continuous calender of all sites and times for merging
  all.sites <- unique(deer.dat$SITE)
  
# create a sequence of date-times from the start to the finish of the survey  
  min.date.time <- as.POSIXct(strptime("2016-06-08 18:00:00", "%Y-%m-%d %H:%M:%S"))
  max.date.time <- as.POSIXct(strptime("2019-01-01 18:00:00", "%Y-%m-%d %H:%M:%S"))  # 2016-06-08 18:00
  # 2018-07-19 15:00
  
# the survey began on Jan 1st at midnight, thus the start date
  seq.date.time <- seq.POSIXt(from = min.date.time, to = max.date.time, by = "hour")
  rep.date.time <- rep(seq.date.time, length(all.sites))
  
# create site reps for new data frame
  rep.sites <- sort(rep(all.sites, length(seq.date.time)))
  
# create data frame for merge
  survey.sites.df <- data.frame(SITE = rep.sites, DATE.HOUR = rep.date.time)
  head(survey.sites.df)
  
# merge site by survey df with aggregated count data
  deer.df <- merge(x = survey.sites.df, y = count.by.hour, by.x = c("SITE","DATE.HOUR"), by.y = c("SITE","DATE.HOUR"), all.x = TRUE, all.y = FALSE)
  deer.df[is.na(deer.df)] <- 0
  head(deer.df, 500)
  tail(deer.df, 500)
  
# prepare backbone for aggregation by hour  
  days.on <- as.numeric(round(max.date.time - min.date.time, 0))
  temp <- deer.df[order(deer.df$DATE.HOUR),]
  rownames(temp) <- 1:nrow(temp)
  temp <- temp[order(temp$SITE, temp$DATE.HOUR),]
  rownames(temp) <- 1:nrow(temp)
  temp$DATE <- as.Date(temp$DATE.HOUR)
  temp <- temp[1:494736,]
  
# add new column to aggregate data by day for modelling population changes
  temp$AGG <- sort(rep(1:days.on, 22*24))
  
# aggregate data based #on AGG column
  temp.agg <- aggregate(COUNT ~ SITE + DATE, data = temp, FUN = "max")
  temp.agg <- temp.agg[order(temp.agg$SITE, temp.agg$DATE),]
  deer.df <- temp.agg
  deer.df$DATE.AGG <- paste(deer.df$DATE, deer.df$AGG, sep=".")
  
# subset count dataframe by month for change in density analysis
  JUN16 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2016")
  JUL16 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2016")
  AUG16 <- subset(deer.df, format.Date(DATE, "%m")=="08"  & format.Date(DATE, "%Y")=="2016")
  OCT16 <- subset(deer.df, format.Date(DATE, "%m")=="10"& format.Date(DATE, "%Y")=="2016" )
  SEP16 <- subset(deer.df, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2016") 
  NOV16 <- subset(deer.df, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2016")
  DEC16 <- subset(deer.df, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2016")
  JAN17 <- subset(deer.df, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2017")
  FEB17 <- subset(deer.df, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2017")
  MAR17 <- subset(deer.df, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2017")
  APR17 <- subset(deer.df, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2017")
  MAY17 <- subset(deer.df, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2017")
  JUN17 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2017")
  JUL17 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2017")
  AUG17 <- subset(deer.df, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2017")
  SEP17 <- subset(deer.df, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2017")
  OCT17 <- subset(deer.df, format.Date(DATE, "%m")=="10" & format.Date(DATE, "%Y")=="2017")
  NOV17 <- subset(deer.df, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2017")
  DEC17 <- subset(deer.df, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2017") 
  JAN18 <- subset(deer.df, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2018")
  FEB18 <- subset(deer.df, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2018")
  MAR18 <- subset(deer.df, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2018")
  APR18 <- subset(deer.df, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2018")
  MAY18 <- subset(deer.df, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2018")
  JUN18 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2018")
  JUL18 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2018")
  AUG18 <- subset(deer.df, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2018")

# rbind into seasons
  summer16 <- rbind(JUN16, JUL16, AUG16)
  fall16 <- rbind(SEP16, OCT16, NOV16)
  winter16 <- rbind(DEC16, JAN17, FEB17)
  spring17 <- rbind(MAR17, APR17, MAY17)
  summer17 <- rbind(JUN17,JUL17, AUG17)
  fall17 <- rbind(SEP17, OCT17, NOV17)
  winter17 <- rbind(DEC17, JAN18, FEB18)
  spring18 <- rbind(MAR18, APR18, MAY18)
  summer18 <- rbind(JUN18,JUL18, AUG18)
  
# create a list of seasons to run functions on simultaneously
  season.list <- list(summer16, fall16, winter16, spring17, summer17, fall17, winter17, spring18, summer18)
  
# create simplified df of counts
  season.list <- lapply(season.list, function(x) x <- x[, c("SITE","DATE","COUNT")])
  
# cast for unmarked  
  season.list <- lapply(season.list, function (x) x <- cast(data = x, formula = SITE~DATE, mean))   #timemeans <- cast(mdata, time~variable, mean) 
  season.list <- lapply(season.list, function (x) x <- as.matrix.cast_df(x[,-1])) 
  
# create site covariate dataframe
  deer.covars <- data.frame(ELE=cam.dat$site.ele, ASPECT = cam.dat$aspect, SLOPE = cam.dat$slope, EDGE=cam.dat$EDGE, ROAD=cam.dat$ROAD)

# create unmarked frames with data 
  season.list <- lapply(season.list, function(x) unmarkedFramePCount(y = x, siteCovs = deer.covars))
  
# model of Royle 2004 based on counts
# general formula is: model <- pcount(~detection_formula ~occupancy_formula, dataframe, K=100, se=TRUE)
# same starts are used from earlier analyses
  model.list <- lapply(season.list, function(x) x <- pcount(~1 ~ASPECT, data = x, K = 150, se = T, mixture = "P",
                                                            starts = c(2.4, 0.005, -6.26),
                                                            control = list(trace = T, REPORT = 1)))
  
 
  
# extract seasonal estimates
  re.list <- lapply(model.list, function(x) ranef(x))
  EBUP<- lapply(re.list, function (x) bup(x, stat = 'mean'))
  CI <- lapply(re.list, function(x) confint(x, level = 0.95))
  
  summer16.EBUP <- unlist(EBUP[[1]])
  summer16.CI <- unlist(CI[[1]])
  summer16.CI <- as.data.frame(summer16.CI)
  summer16.EST <- data.frame(Estimate = sum(summer16.EBUP), Lower = sum(summer16.CI$`2.5%`), Upper = sum(summer16.CI$`97.5%`))
  
  fall16.EBUP <- unlist(EBUP[[2]])
  fall16.CI <- unlist(CI[[2]])
  fall16.CI <- as.data.frame(fall16.CI)
  fall16.EST <- data.frame(Estimate = sum(fall16.EBUP), Lower = sum(fall16.CI$`2.5%`), Upper = sum(fall16.CI$`97.5%`))
  
  winter16.EBUP <- unlist(EBUP[[3]])
  winter16.CI <- unlist(CI[[3]])
  winter16.CI <- as.data.frame(winter16.CI)
  winter16.EST <- data.frame(Estimate = sum(winter16.EBUP), Lower = sum(winter16.CI$`2.5%`), Upper = sum(winter16.CI$`97.5%`))
  
  spring17.EBUP <- unlist(EBUP[[4]])
  spring17.CI <- unlist(CI[[4]])
  spring17.CI <- as.data.frame(spring17.CI)
  spring17.EST <- data.frame(Estimate = sum(spring17.EBUP), Lower = sum(spring17.CI$`2.5%`), Upper = sum(spring17.CI$`97.5%`))
  
  summer17.EBUP <- unlist(EBUP[[5]])
  summer17.CI <- unlist(CI[[5]])
  summer17.CI <- as.data.frame(summer17.CI)
  summer17.EST <- data.frame(Estimate = sum(summer17.EBUP), Lower = sum(summer17.CI$`2.5%`), Upper = sum(summer17.CI$`97.5%`))
  
  fall17.EBUP <- unlist(EBUP[[6]])
  fall17.CI <- unlist(CI[[6]])
  fall17.CI <- as.data.frame(fall17.CI)
  fall17.EST <- data.frame(Estimate = sum(fall17.EBUP), Lower = sum(fall17.CI$`2.5%`), Upper = sum(fall17.CI$`97.5%`))
  
  winter18.EBUP <- unlist(EBUP[[7]])
  winter18.CI <- unlist(CI[[7]])
  winter18.CI <- as.data.frame(winter18.CI)
  winter18.EST <- data.frame(Estimate = sum(winter18.EBUP), Lower = sum(winter18.CI$`2.5%`), Upper = sum(winter18.CI$`97.5%`))
  
  spring18.EBUP <- unlist(EBUP[[8]])
  spring18.CI <- unlist(CI[[8]])
  spring18.CI <- as.data.frame(spring18.CI)
  spring18.EST <- data.frame(Estimate = sum(spring18.EBUP), Lower = sum(spring18.CI$`2.5%`), Upper = sum(spring18.CI$`97.5%`))
  
  summer18.EBUP <- unlist(EBUP[[9]])
  summer18.CI <- unlist(CI[[9]])
  summer18.CI <- as.data.frame(summer18.CI)
  summer18.EST <- data.frame(Estimate = sum(summer18.EBUP), Lower = sum(summer18.CI$`2.5%`), Upper = sum(summer18.CI$`97.5%`))
  

# bind estimates togehter
  month.est <- rbind(summer16.EST, fall16.EST, winter16.EST, spring17.EST, summer17.EST, fall17.EST, winter18.EST, spring18.EST, summer18.EST)
  month.est <- month.est/10.24  
  names(month.est) <- c("Estimate", "Lower", "Upper")
  month.est$SEASON <- c('summer16', 'fall16', 'winter16', 'spring17', 'summer17', 'fall17', 'winter17', 'spring18', 'summer18')
  
# read-in ndvi season data
  ndvi.dat <- read.csv("ndvi_full_period.csv")
  
# convert date to date
  ndvi.dat$DATE <- as.Date(ndvi.dat$DATE)
  
# subset ndvi by season
  JUN16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2016")
  JUL16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2016")
  AUG16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08"  & format.Date(DATE, "%Y")=="2016")
  OCT16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="10"& format.Date(DATE, "%Y")=="2016")
  SEP16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2016")
  NOV16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2016")
  DEC16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2016")
  JAN17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2017")
  FEB17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2017")
  MAR17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2017")
  APR17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2017")
  MAY17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2017")
  JUN17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2017")
  JUL17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2017")
  AUG17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2017")
  SEP17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2017")
  OCT17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="10" & format.Date(DATE, "%Y")=="2017")
  NOV17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2017")
  DEC17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2017") 
  JAN18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2018")
  FEB18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2018")
  MAR18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2018")
  APR18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2018")
  MAY18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2018")
  JUN18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2018")
  JUL18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2018")
  AUG18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2018")
  
# rbind into seasons
  summer16.ndvi <- rbind(JUN16, JUL16, AUG16)
  fall16.ndvi <- rbind(SEP16, OCT16, NOV16)
  winter16.ndvi <- rbind(DEC16, JAN17, FEB17)
  spring17.ndvi <- rbind(MAR17, APR17, MAY17)
  summer17.ndvi <- rbind(JUN17,JUL17, AUG17)
  fall17.ndvi <- rbind(SEP17, OCT17, NOV17)
  winter17.ndvi <- rbind(DEC17, JAN18, FEB18)
  spring18.ndvi <- rbind(MAR18, APR18, MAY18)
  summer18.ndvi <- rbind(JUN18,JUL18, AUG18)
 
  
# add factor for season
  summer16.ndvi$SEASON <- rep(1, nrow(summer16.ndvi))
  sum16mean <- mean(summer16.ndvi$NDVI)
  bstrap <- NULL
  for (i in 1:10000){
    bstrap <- c(bstrap, mean(sample(summer16.ndvi$NDVI, nrow(summer16.ndvi),replace=T)))}
  ndvisum16 <- data.frame(Mean = sum16mean, Lower = quantile(bstrap,.025), Upper =quantile(bstrap, .975))
  
  
  fall16.ndvi$SEASON <- rep(2, nrow(fall16.ndvi))
  fstrap <- NULL
  for (i in 1:10000){
    fstrap <- c(fstrap, mean(sample(fall16.ndvi$NDVI, nrow(fall16.ndvi),replace=T)))}
  ndvifall16 <- data.frame(Mean = mean(fall16.ndvi$NDVI), Lower = quantile(fstrap,.025), Upper =quantile(fstrap, .975))
  
  
  winter16.ndvi$SEASON <- rep(3, nrow(winter16.ndvi))
  wstrap <- NULL
  for (i in 1:10000){
    wstrap <- c(wstrap, mean(sample(winter16.ndvi$NDVI, nrow(winter16.ndvi),replace=T)))}
  ndviwinter16 <- data.frame(Mean = mean(winter16.ndvi$NDVI), Lower = quantile(wstrap,.025), Upper =quantile(wstrap, .975))
  
  
  spring17.ndvi$SEASON <-rep(4, nrow(spring17.ndvi))
  sprstrap <- NULL
  for (i in 1:10000){
    sprstrap <- c(sprstrap, mean(sample(spring17.ndvi$NDVI, nrow(spring17.ndvi),replace=T)))}
  ndvispring17 <- data.frame(Mean = mean(spring17.ndvi$NDVI), Lower = quantile(sprstrap,.025), Upper =quantile(sprstrap, .975))
  
  
  summer17.ndvi$SEASON <- rep(5, nrow(summer17.ndvi))
  sumstrap <- NULL
  for (i in 1:10000){
    sumstrap <- c(sumstrap, mean(sample(summer17.ndvi$NDVI, nrow(summer17.ndvi),replace=T)))}
  ndvisummer17 <- data.frame(Mean = mean(summer17.ndvi$NDVI), Lower = quantile(sumstrap,.025), Upper =quantile(sumstrap, .975))
  
  
  fall17.ndvi$SEASON <- rep(6, nrow(fall17.ndvi))
  fastrap <- NULL
  for (i in 1:10000){
    fastrap <- c(fastrap, mean(sample(fall17.ndvi$NDVI, nrow(fall17.ndvi),replace=T)))}
  ndvifall17 <- data.frame(Mean = mean(fall17.ndvi$NDVI), Lower = quantile(fastrap,.025), Upper =quantile(fastrap, .975))
  
  winter17.ndvi$SEASON <- rep(7, nrow(winter17.ndvi))
  wstrap <- NULL
  for (i in 1:10000){
    wstrap <- c(wstrap, mean(sample(winter17.ndvi$NDVI, nrow(winter17.ndvi),replace=T)))}
  ndviwinter17 <- data.frame(Mean = mean(winter17.ndvi$NDVI), Lower = quantile(wstrap,.025), Upper =quantile(wstrap, .975))
  
  spring18.ndvi$SEASON <-rep(8, nrow(spring18.ndvi))
  sprstrap <- NULL
  for (i in 1:10000){
    sprstrap <- c(sprstrap, mean(sample(spring18.ndvi$NDVI, nrow(spring18.ndvi),replace=T)))}
  ndvispring18 <- data.frame(Mean = mean(spring18.ndvi$NDVI), Lower = quantile(sprstrap,.025), Upper =quantile(sprstrap, .975))
  
  
  summer18.ndvi$SEASON <- rep(9, nrow(summer18.ndvi))
  sumstrap <- NULL
  for (i in 1:10000){
    sumstrap <- c(sumstrap, mean(sample(summer18.ndvi$NDVI, nrow(summer18.ndvi),replace=T)))}
  ndvisummer18 <- data.frame(Mean = mean(summer18.ndvi$NDVI), Lower = quantile(sumstrap,.025), Upper =quantile(sumstrap, .975))
  
# rebind
  ndvi.dat <- rbind(summer16.ndvi, fall16.ndvi, winter16.ndvi, spring17.ndvi, summer17.ndvi, fall17.ndvi, winter17.ndvi, spring18.ndvi, summer18.ndvi)
  ndvi.plot <- rbind(ndvisum16,ndvifall16, ndviwinter16, ndvispring17, ndvisummer17, ndvifall17, ndviwinter17, ndvispring18, ndvisummer18)

# plot the population changes with NDVI
  y <- month.est$Estimate
  ndvi <- ndvi.plot$Mean / 10000
  ndvi <- ndvi*40
  ndvi <- matrix(ndvi,1,10)
  
  
  
  df.bar <- barplot(ndvi, xlab = 'Season', ylab = "Deer Density (Deer / Sq.Km)", ylim = c(0, 40), col = 'light green', space = 0, width = 1, cex.lab = 1.5, cex.axis = 1.5)
  
  
  axis(side = 4, cex.axis = 1.5, col.axis = 'forest green', labels = c(0, 0.25, 0.5, 0.75, 1), at = c(0, 10 , 20 , 30, 40))
  mtext("NDVI", side = 4, cex = 1.5, line = 3)
  
  
  
  points(y ~ df.bar,
         xlab = "Season",
         ylab = "Density (deer / sq.km)",
         cex.lab = 1.5,
         cex.axis = 1.5,
         pch = 19,
         ylim = c(0,40),
         xaxt = 'n',
         cex = 1.5
  )
  
  
  axis(side = 1, at = df.bar, c('Summer 2016', 'Fall 2016', 'Winter 2016', 'Spring 2017', 'Summer 2017', 'Fall 2017', 'Winter 2017', 'Spring 2018', 'Summer 2018', 'Fall 2018'), cex = 1.5)
  arrows(x0 = df.bar, x1 = df.bar, y0 = month.est$Estimate, y1 = month.est$Upper, angle = 90, length = 0.05)
  arrows(x0 = df.bar, x1 = df.bar, y0 = month.est$Estimate, y1 = month.est$Lower, angle = 90, length = 0.05)
  
  arrows(x0 = 7, x1 = 7, y0 = 31, y1 = 40, angle = 90, length = 0.05)  
  arrows(x0 = 7, x1 = 7, y0 = 31, y1 = 23.7, angle = 90, length = 0.05) 
  
  lines(df.bar, y, col='black', lwd=2)
 
  points(31 ~ 7,
  pch = 19,
  cex = 1.5,
  col = 'dark grey')
  
  
  