require(CircStats)
require(adehabitatHR)

anglefun <- function(xx,yy,bearing=TRUE,as.deg=FALSE){
  ## calculates the compass bearing of the line between two points
  ## xx and yy are the differences in x and y coordinates between two points
  
  ## Options:
  ## bearing = FALSE returns +/- pi instead of 0:2*pi
  ## as.deg = TRUE returns degrees instead of radians
  c = 1
  if (as.deg){
    c = 180/pi
  }
  
  b<-sign(xx)
  b[b==0]<-1  #corrects for the fact that sign(0) == 0
  tempangle = b*(yy<0)*pi+atan(xx/yy)
  if(bearing){
    #return a compass bearing 0 to 2pi
    #if bearing==FALSE then a heading (+/- pi) is returned
    tempangle[tempangle<0]<-tempangle[tempangle<0]+2*pi
  }
  return(tempangle*c)
}

bearing.ta <- function(locs,as.deg=FALSE){
  ## calculates the bearing and length of the two lines
  ##    formed by three points
  ## the turning angle from the first bearing to the
  ##    second bearing is also calculated
  ## locations are assumed to be in (X,Y) format.
  ## Options:
  ## as.deg = TRUE returns degrees instead of radians
  if (length(locs[1:2]) != 2 | length(locs[3:4]) != 2 | length(locs[5:6]) !=2){
    print("Locations must consist of either three vectors, length == 2,
or three two-column dataframes")
    return(NaN)
  }
  
  locdiff1<-locs[3:4]-locs[1:2]
  locdiff2<-locs[5:6]-locs[3:4]
  bearing1<-anglefun(locdiff1[1],locdiff1[2],bearing=F) %>% .[,1]
  bearing2<-anglefun(locdiff2[1],locdiff2[2],bearing=F) %>% .[,1]
  
  ta=(bearing2-bearing1)
  
  ta1=head(ta, -1)
  ta1[ which(ta1 < -pi) ] <- ta1[ which(ta1 < -pi) ] + 2*pi
  ta1[ which(ta1 > pi) ] <- ta1[ which(ta1 > pi) ] - 2*pi
  return(ta1)
}

bearing.sl <- function(locs,as.deg=FALSE){
  ## calculates the bearing and length of the two lines
  ##    formed by three points
  ## the turning angle from the first bearing to the
  ##    second bearing is also calculated
  ## locations are assumed to be in (X,Y) format.
  ## Options:
  ## as.deg = TRUE returns degrees instead of radians
  if (length(locs[1:2]) != 2 | length(locs[3:4]) != 2 | length(locs[5:6]) !=2){
    print("Locations must consist of either three vectors, length == 2,
or three two-column dataframes")
    return(NaN)
  }
  
  locdiff1<-locs[3:4]-locs[1:2]
  locdiff2<-locs[5:6]-locs[3:4]
  
  if(is.data.frame(locdiff1)){
    dist1<-sqrt(rowSums(locdiff1^2))
    dist2<-sqrt(rowSums(locdiff2^2))
  }else{
    dist1<-sqrt(sum(locdiff1^2))
    dist2<-sqrt(sum(locdiff2^2))
  }
  
  return(dist2)
}

format.ssf=function(datfr){
  datfr = datfr %>% arrange(Time) %>% dplyr::distinct(Time, .keep_all = TRUE)
  ## Create a new dataframe with a starting and ending time and location for each step
  stps=with(datfr, data.frame(# Duplicate time column at a lag of 1 as a new column
    Time_1=c(Time), Time_2=c(Time[-1], NA), Time_3=c(Time[c(-1,-2)], NA, NA),
    # Duplicate lat and long columns at a lag of 1 as new columns
    Lat_1=c(Lat), Lon_1=c(Lon), Lat_2=c(Lat[-1], NA), Lon_2=c(Lon[-1],NA),
    # Duplicate lat and long columns at a lag of 2 as new columns
    Lat_3=c(Lat[-c(1:2)], NA,NA), Lon_3=c(datfr$Lon[-c(1:2)],NA,NA),
    # Create two SL and TA columns populated with NA's
    SL=rep(NA, nrow(datfr)), TA=rep(NA, nrow(datfr)),
    # Carry through calf status
    CalfStatus=c(CalfStatus[-c(1:2)], NA,NA),
    # Create a column indicating these are the "true" steps and not generated ones
    Type = rep(TRUE, nrow(datfr))))
  
  stps$TimeDiff = stps$Time_3 - stps$Time_2
  
  ## Calculates the TA's and SL's for each set of points 
  stps$TA=c(bearing.ta(stps[,c("Lat_1", "Lon_1","Lat_2", "Lon_2","Lat_3", "Lon_3")], as.deg=FALSE), NA)
  # Adds a SL and TA column to the dataframe
  stps$SL=c(bearing.sl(stps[,c("Lat_1", "Lon_1","Lat_2", "Lon_2","Lat_3", "Lon_3")], as.deg=FALSE))
  
  # stps$f_ch = ifelse(stps$Lat_1 == stps$Lat_2 & stps$Lon_1 == stps$Lon_2, FALSE, TRUE)
  # stps$s_ch = ifelse(stps$Lat_2 == stps$Lat_3 & stps$Lon_2 == stps$Lon_3, FALSE, TRUE)
  
  ## Removes last two lines of dataframe as they don't have TA's of SL's
  stps=slice(subset(stps, select = -c(Lat_1, Lon_1, Time_1)), 1:(n()-2)) %>% 
    dplyr::rename(Lat_1 = Lat_2, Lon_1 = Lon_2, Lat_2 = Lat_3, Lon_2 = Lon_3,
                  Time_1 = Time_2, Time_2 = Time_3, CalfStatus = CalfStatus)
  return(stps)
}


# gen.stps=function(stps, stepnum){
#   ## Set the number of random steps to be generated
#   nstps = stepnum
#   
#   
#   ## Fit data to a gamma distribution to get parameter estimates
#   # Subset data by those calculated without skipped fixes
#   dat.sl=(dplyr::filter(stps, TimeDiff==median(TimeDiff))$SL)/1000
#   
#   # Create the likelihood function for a gamma dist
#   gamma_loglik=function(parm){
#     shape=parm[1]
#     rate=parm[2]
#     loglik=sum(log(dgamma(dat.sl, shape, rate)))
#     return(-loglik)}
#   # Create object with parameter estimates
#   suppressWarnings(sl.est<-nlminb(parm<-c(1,1), gamma_loglik, hessian=TRUE))
#   
#   
#   
#   ## Fit data to a von-mises distribution to get parameter estimates
#   # Omit data with skipped fixes, difference fix schedules, etc
#   dat.ta=(subset(stps, TimeDiff==median(TimeDiff)) %>%
#             # Omit data with step lengths lower than 10m as these alter turning angle dists
#             dplyr::filter(SL >= 10))$TA
#   # Create the likelihood function for a von mises dist
#   vonmis_loglik=function(parm){
#     mu=parm[1]
#     kappa=parm[2]
#     loglik=sum(log(dvm(dat.ta, mu, kappa)))
#     return(-loglik)}
#   # Create object with parameter estimates
#   suppressWarnings(ta.est<-nlminb(parm<-c(1,2), vonmis_loglik, hessian=TRUE))
#   
#   
#   ## Replacing super low parameter estimates with minimum threshold for dvm function
#   ta.est$par[2] = ifelse(ta.est$par[2] < 0.0000001, 0.000001, ta.est$par[2])
#   
#   ## Create column indicating these are comparison steps
#   stps1 = stps %>% subset(select = -c(SL, TA, Type)) %>%
#     ## Remove the SL and TA data from real steps (to be replaced with generated SL and TA)
#     add_column(Type = rep(FALSE, nrow(.))) %>%
#     ## Repeat each line for the number of steps to be generated
#     add_column(Count = rep(nstps, nrow(.))) %>%
#     uncount(Count)
#   
#   
#   ## Generate random step lengths from a gamma distribution with parameters esitmated above
#   stps1$SL=sample
#   ## Generate random turning angles from a vonmises dist. with parameters estimated above
#   stps1$TA=(rvm(nrow(stps1), ta.est$par[1], ta.est$par[2])) %%
#     # Convert these TA's from 0 - 2pi to pi to -pi
#     (2*pi) %>%
#     ifelse(. > pi, . - (2* pi), .)
#   ## Calculate the longitude values of the new steps
#   stps1$Lon_2<-stps1$Lon_1 + stps1$SL*cos(stps1$TA)
#   ## Calculate the latitude values of the new steps
#   stps1$Lat_2<-stps1$Lat_1 + stps1$SL*sin(stps1$TA)
#   
#   
#   ## Bind dataframes of true and generated steps together
#   steps=rbind(stps,stps1 %>% dplyr::select(-Type, -TimeDiff, everything())) %>%
#     # Arrange new dataframe by timestep
#     arrange(Time_1) %>%
#     # Adding a column to group individual steps together with their replicates
#     add_column(Step = paste(rep('step_', nrow(stps1)), rep(1:nrow(stps), 1, each= (nstps+1)), sep = ""))
#   
#   return(steps)
# }



gen.stps=function(stps, stepnum){
  ## Set the number of random steps to be generated
  nstps = stepnum
  
  
  ## Fit data to a gamma distribution to get parameter estimates
  # Subset data by those calculated without skipped fixes
  dat.sl=(dplyr::filter(stps, TimeDiff==median(TimeDiff))$SL)*1000

  # Create the likelihood function for a gamma dist
  gamma_loglik=function(parm){
    shape=parm[1]
    rate=parm[2]
    loglik=sum(log(dgamma(dat.sl, shape, rate)))
    return(-loglik)}
  # Create object with parameter estimates
  suppressWarnings(sl.est<-nlminb(parm<-c(1,1), gamma_loglik, hessian=TRUE))


  
  ## Fit data to a von-mises distribution to get parameter estimates
  # Omit data with skipped fixes, difference fix schedules, etc
  dat.ta=(subset(stps, TimeDiff==median(TimeDiff)) %>%
            # Omit data with step lengths lower than 10m as these alter turning angle dists
            dplyr::filter(SL >= 0.00001))$TA
  # Create the likelihood function for a von mises dist
  vonmis_loglik=function(parm){
    mu=parm[1]
    kappa=parm[2]
    loglik=sum(log(dvm(dat.ta, mu, kappa)))
    return(-loglik)}
  # Create object with parameter estimates
  suppressWarnings(ta.est<-nlminb(parm<-c(1,2), vonmis_loglik, hessian=TRUE))

  
  ## Replacing super low parameter estimates with minimum threshold for dvm function
  ta.est$par[2] = ifelse(ta.est$par[2] < 0.0000001, 0.000001, ta.est$par[2])
  
  ## Create column indicating these are comparison steps
  stps1 = stps %>% subset(select = -c(SL, TA, Type)) %>%
    ## Remove the SL and TA data from real steps (to be replaced with generated SL and TA)
    add_column(Type = rep(FALSE, nrow(.))) %>%
    ## Repeat each line for the number of steps to be generated
    add_column(Count = rep(nstps, nrow(.))) %>%
    uncount(Count)
  
  
  ## Generate random step lengths from a gamma distribution with parameters esitmated above
  stps1$SL=rgamma(nrow(stps1), sl.est$par[1], sl.est$par[2])/1000
  ## Generate random turning angles from a vonmises dist. with parameters estimated above
  stps1$TA=(rvm(nrow(stps1), ta.est$par[1], ta.est$par[2])) %%
    # Convert these TA's from 0 - 2pi to pi to -pi
    (2*pi) %>%
    ifelse(. > pi, . - (2* pi), .)
  ## Calculate the longitude values of the new steps
  stps1$Lon_2<-stps1$Lon_1 + stps1$SL*cos(stps1$TA)
  ## Calculate the latitude values of the new steps
  stps1$Lat_2<-stps1$Lat_1 + stps1$SL*sin(stps1$TA)
  
  
  ## Bind dataframes of true and generated steps together
  steps=rbind(stps,stps1 %>% dplyr::select(-Type, -TimeDiff, everything())) %>%
    # Arrange new dataframe by timestep
    arrange(Time_1) %>%
    # Adding a column to group individual steps together with their replicates
    add_column(Step = paste(rep('step_', nrow(stps1)), rep(1:nrow(stps), 1, each= (nstps+1)), sep = ""))
  
  return(steps)
}


## Function calculating dispersal from 5th point 
ann_disp = function(df){
  # Isolating fifth point and getting standard lat/long
  init_pt = df %>% dplyr::filter(row_number() == 1) %>% dplyr::select(Long_std, Lat_std)
  # Isolating all other standard lat/longs
  all_pts = df %>% dplyr::select(Long_std, Lat_std) 
  # Returning vector of displacement distances in meters
  c(distm(all_pts, init_pt, fun = distHaversine)/1000)
}


## Function to fix intercept variance in models before running them
var.fix=function(TMBobj){
  ## Fixing the standard deviation for the random ID component to a large number
  TMBobj$parameters$theta[1] = log(1e3)
  ## Determining the number of parameters to allow estimation in the list below
  npar = length(TMBobj$parameters$theta)
  ## Allow all other variances to be estimated except the first one 
  TMBobj$mapArg = list(theta=factor(c(NA, rep(1, (npar-1)))))
  return(TMBobj)
}



near_wlf_ssf = function(datfr_pry, datfr_prd, type){
  if(type == "end"){
  ref <- datfr_pry$Time_2
  
  dists <- map_dbl(seq_along(ref), \(i) {
    tme  <- ref[i]
    foc  <- datfr_pry[i, c("Lon_2", "Lat_2")]
    datfr_prd %>% 
      mutate(diff_h = abs(as.numeric(difftime(Time, tme, "hours")))) %>% 
      dplyr::filter(diff_h <= 5) %>% 
      group_by(AnimalId) %>% 
      slice_min(diff_h, with_ties = FALSE) %>% 
      mutate(dist_km = geosphere::distHaversine(cbind(Lon, Lat), foc)/1000) %>% ungroup() %>%
      #summarise(med_km = median(dist_km)) %>%   # median here
      pull(dist_km) %>% min()
  })
  }
  
  
  if(type == "start"){
    ref <- datfr_pry %>% filter(Type == TRUE) %>% .$Time_1
  
  dists <- map_dbl(seq_along(ref), \(i) {
    tme  <- ref[i]
    foc  <- datfr_pry %>% filter(Type == TRUE) %>% .[i, c("Lon_1", "Lat_1")]
    datfr_prd %>% 
      mutate(diff_h = abs(as.numeric(difftime(Time, tme, "hours")))) %>% 
      dplyr::filter(diff_h <= 5) %>% 
      group_by(AnimalId) %>% 
      slice_min(diff_h, with_ties = FALSE) %>% 
      mutate(dist_km = geosphere::distHaversine(cbind(Lon, Lat), foc)/1000) %>% ungroup() %>%
      #summarise(med_km = median(dist_km)) %>%   # median here
      pull(dist_km) %>% min()
  })
}
  if(type == "end"){
  datfr_pry = datfr_pry %>% mutate(dst_min_end = dists)
  }
  if(type == "start"){
  datfr_pry = datfr_pry %>% mutate(dst_min_strt = rep(dists, each = 16))
  }
  
  return(datfr_pry)
}

near_wlf_rsf = function(datfr_pry, datfr_prd){

    ref <- datfr_pry$Time
    
    
    dists <- map_dbl(seq_along(ref), \(i) {
      tme  <- ref[i]
      foc  <- datfr_pry[i, c("Lon", "Lat")]
      datfr_prd %>% 
        mutate(diff_h = abs(as.numeric(difftime(Time, tme, "hours")))) %>% 
        dplyr::filter(diff_h <= 5) %>% 
        group_by(AnimalId) %>% 
        slice_min(diff_h, with_ties = FALSE) %>% 
        mutate(dist_km = geosphere::distHaversine(cbind(Lon, Lat), foc)/1000) %>% ungroup() %>%
        #summarise(med_km = median(dist_km)) %>%   # median here
        pull(dist_km) %>% min()
    })

  
  

    datfr_pry = datfr_pry %>% mutate(dst_min = dists)
  
  return(datfr_pry)
}








