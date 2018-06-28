# Modelling the historical deforestation on the island of Mauritius, 
# using a cellular automata approach.
 
# This model simulates the historical deforestation on the island of Mauritius using six available 
# historical maps and several other suitability maps.
 
# Created by: Matthijs Hinkamp, University of Amsterdam

# Before running, make sure to change the directory and to put all of the .txt files in 
# the ready data map directly into that directory.

rm(list=ls())

## Input
###########################################
start_year = 1638 # options: 1638, 1685, 1773, 1835, 1872, 1935 and 1997 (for now only 1638)
end_year = 1997 # options: 1686, 1774, 1836, 1873, 1936 and 1997

# Choose deforestation rate calculation method
# The default is a linear interpolation between the historical maps.
# If improved rate is on (improved_rate==1) a rate will be used which is partially backed up by historical sources.
improved_rate=1

# Check rate?
# This function checks amount of simulated forest cells and the amount of forest cells for the known maps. 
# Use this when changing the model to check proper functioning. Acceptable values are -20<y<20.
check_rate=1

# Print accuracy?
# This fuction will print the accuracy of the final results. This works when calubration is on and off.
print_accuracy=1

# Choose the suitability type
# Selecting historical suitability changes the suitability per period based on the known locations of settlements.
# If historical suitability is off, a static suitability is assumed based on environmental factors.
historical_suitability=1

# Calibration?
# Calibration goes through a range of weights for one factor at a time and returns an array (and graph) with the accuracy for each value
# Calibration requires historical suitability and improved deforestation rate.
# If calibration is on, make sure the visualisations are turned off.
# Calibration runs a number of calibration runs with different weights
calibration=1
# Initial weights
inw=50
# Number of calibration runs (the first run is with the initial values)
calibration_runs=101
# Calibration steps
# How much should the weights shift per iteration?
calibration_steps=1
# Set (initial) weights manually
# Dutch occupation
w1_sett_dutch=inw
w1_dem=0.0001
w1_slope=0.0001
w1_agr_suitab=0.0001

# French occupation
w2_sett_french=53
w2_riv_french=72
w2_dem=0.0001
w2_slope=0.0001
w2_agr_suitab=3

# English occupation (early)
w3_coast_dist=0.0001
w3_dem=75
w3_slope=7
w3_agr_suitab=0.0001

# English occupation (mid)
w4_coast_dist=7
w4_dem=57
w4_slope=0.0001
w4_agr_suitab=58

# English occupation (late)
w5_coast_dist=0.0001
w5_dem=53
w5_slope=22
w5_agr_suitab=73

# Indpendence
w6_coast_dist=0.0001
w6_dem=inw
w6_slope=1
w6_agr_suitab=0.0001

# Which weight has to be calibrated?
cal_w1_sett_dutch=0
cal_w1_dem=0
cal_w1_slope=0
cal_w1_agr_suitab=0

cal_w2_sett_french=0
cal_w2_dem=0
cal_w2_slope=0
cal_w2_agr_suitab=0
cal_w2_riv_french=0

cal_w3_coast_dist=0
cal_w3_dem=0
cal_w3_slope=0
cal_w3_agr_suitab=0

cal_w4_coast_dist=0
cal_w4_dem=0
cal_w4_slope=0
cal_w4_agr_suitab=0

cal_w5_coast_dist=0
cal_w5_dem=0
cal_w5_slope=0
cal_w5_agr_suitab=0

cal_w6_coast_dist=1
cal_w6_dem=0
cal_w6_slope=0
cal_w6_agr_suitab=0

# Choose visualisation (0=off and 1=on)
graph_deforestation_rate=0
graph_forest_extent=0

graph_regrowth=0

images_forest_extent=0

dynamic_visualisation=0

## Initialisation
###########################################

v=0
c=1
cal_st <- c(-calibration_steps,calibration_steps)
forest_record=0
rate_adj_record=0
rotate <- function(x) t(apply(x, 2, rev))
mycols <- colors()[c(123,649,139)]
regrowth <- matrix(nrow = 602,ncol = 524)
regrowth[1:315448]=0
regrowth_record=0

acc_1685=0
acc_1773=0
acc_1835=0
acc_1872=0
acc_1935=0
acc_1997=0
acc_engind=0

# Directory
getwd()
setwd(dir <- "C:/users/Matthijs/r_directory")

# Loading input data
forest_1638 <- data.matrix(read.table("forest_1638.txt"),rownames.force = NA)
forest_1685 <- data.matrix(read.table("forest_1685.txt"),rownames.force = NA)
forest_1773 <- data.matrix(read.table("forest_1773.txt"),rownames.force = NA)
forest_1835 <- data.matrix(read.table("forest_1835.txt"),rownames.force = NA)
forest_1872 <- data.matrix(read.table("forest_1872.txt"),rownames.force = NA)
forest_1935 <- data.matrix(read.table("forest_1935.txt"),rownames.force = NA)
forest_1997 <- data.matrix(read.table("forest_1997.txt"),rownames.force = NA)

forest <- forest_1638

agr_suitab_h <- data.matrix(read.table("agr_suitab_h.txt"),rownames.force = NA)
agr_suitab_l <- data.matrix(read.table("agr_suitab_l.txt"),rownames.force = NA)

dem <- data.matrix(read.table("dem.txt"),rownames.force = NA)

dist_coast <- data.matrix(read.table("dist_coast.txt"),rownames.force = NA)

slope_angle <- data.matrix(read.table("slope_angle.txt"),rownames.force = NA)

bat <- data.matrix(read.table("bat.txt"),rownames.force = NA)
br <- data.matrix(read.table("br.txt"),rownames.force = NA)
flacq <- data.matrix(read.table("flacq.txt"),rownames.force = NA)
gp <- data.matrix(read.table("gp.txt"),rownames.force = NA)
sett_1772 <- data.matrix(read.table("sett_1772.txt"),rownames.force = NA)
rivers_1772 <- data.matrix(read.table("rivers2_1772.txt"),rownames.force = NA)
sett_1638 <- data.matrix(read.table("sett_1638.txt"),rownames.force = NA)
ted <- data.matrix(read.table("ted.txt"),rownames.force = NA)
map <- data.matrix(read.table("mp.txt"),rownames.force = NA)
pl <- data.matrix(read.table("pl.txt"),rownames.force = NA)

rows <- nrow(dem)
columns <- ncol(dem)
x <- 1:ncol(dem)
y <- 1:nrow(dem)
cells_island=length(which(forest_1638==2))

# calculating static and initial suitability
suitab_agr_l <- ifelse(agr_suitab_l>0,1-agr_suitab_l/max(agr_suitab_l),agr_suitab_l)
suitab_agr_h <- ifelse(agr_suitab_h>0,1-agr_suitab_h/max(agr_suitab_h),agr_suitab_h)

suitab_dem <- ifelse(dem>0,1-dem/max(dem),dem)

suitab_slope <- ifelse(slope_angle>0,1-slope_angle/max(slope_angle),slope_angle)

suitab_coast_dist <-ifelse(dist_coast>0,1-dist_coast/max(dist_coast),dist_coast)

suitab_bat <- ifelse(bat>0,1-bat/max(bat),bat)
suitab_br <- ifelse(br>0,1-br/max(br),br)
suitab_flacq <- ifelse(flacq>0,1-flacq/max(flacq),flacq)
suitab_gp <- ifelse(gp>0,1-gp/max(gp),gp)
suitab_sett_1638 <- ifelse(sett_1638>0,1-sett_1638/max(sett_1638),sett_1638)
suitab_sett_1772 <- ifelse(sett_1772>0,1-sett_1772/max(sett_1772),sett_1772)
suitab_rivers_1772 <- ifelse(rivers_1772>0,1-rivers_1772/max(rivers_1772),rivers_1772)
suitab_ted <- ifelse(ted>0,1-ted/max(ted),ted)
suitab_pl <- ifelse(pl>0,1-pl/max(pl),pl)
suitab_map <- ifelse(map>0,1-map/max(map),map)

if (historical_suitability==0){
  suitab_total <- (suitab_agr_l+suitab_coast_dist+suitab_dem+suitab_slope)/4
}

## Transition rate (deforestation rate)
###########################################

# Determining the forest extent for the available years
forest_ext_1638 <- length(which(forest_1638==2))
forest_ext_1685 <- length(which(forest_1685==2))
forest_ext_1773 <- length(which(forest_1773==2))
forest_ext_1835 <- length(which(forest_1835==2))
forest_ext_1872 <- length(which(forest_1872==2))
forest_ext_1935 <- length(which(forest_1935==2))
forest_ext_1997 <- length(which(forest_1997==2))

# Linear transition rate (for now, might change later)
rate_1638_1685 <- rep((forest_ext_1638-forest_ext_1685)/(1685-1639),1685-(1638+1))
rate_1685_1773 <- rep((forest_ext_1685-forest_ext_1773)/(1773-1685),1773-1685)
rate_1773_1835 <- rep((forest_ext_1773-forest_ext_1835)/(1835-1773),1835-1773)
rate_1835_1872 <- rep((forest_ext_1835-forest_ext_1872)/(1872-1835),1872-1835)
rate_1872_1935 <- rep((forest_ext_1872-forest_ext_1935)/(1935-1872),1935-1872)
rate_1935_1997 <- rep((forest_ext_1935-forest_ext_1997)/(1997-1935),1997-1935)
rate <- c(rate_1638_1685,rate_1685_1773,rate_1773_1835,rate_1835_1872,rate_1872_1935,rate_1935_1997)

# Creating a discrete transition rate
rate_floor <- floor(rate)
rate_decimals <- (rate-trunc(rate))*10

rounding_random <- runif(1997-(1638+1),0.0,10.0)

for (s in 1:length(rate)){
  rate[s] <- ifelse (rounding_random[s]<rate_decimals[s],(rate_floor[s]+1),rate_floor[s])
}

rate <- append(0,rate)
rate <- append(rate,0)

## The model
###################################################

if (calibration==0){
  # Improved rate (more historical)
  if (improved_rate==1){
    rate2 <- matrix(nrow = 3,ncol=360)
    rate2[1,1:360] = 1638:1997
    rate2[2,1:360] = 1:360
    # 1638
    rate2[3,1]=0
    # 1639-1645
    rate2[3,2:8]=175
    # 1646-1650
    rate2[3,9:13]=220
    # 1651-1655
    rate2[3,14:18]=80
    # 1656-1658
    rate2[3,19:21]=118
    # 1659-1663
    rate2[3,22:26]=0
    # 1664-1668
    rate2[3,27:31]=170
    # 1669-1685
    rate2[3,32:47]=145
    # 1686-1705
    rate2[3,48:69]=100
    # 1706-1721
    rate2[3,70:84]=0
    # 1722-1766
    rate2[3,85:129]=479
    # 1766-1790
    rate2[3,130:153]=697.666
    # 1791-1845
    rate2[3,154:208]=905.9
    # 1846-1870
    rate2[3,209:233]=1647.5
    # 1871-1890
    rate2[3,234:253]=1011
    # 1891-1920
    rate2[3,254:283]=281
    # 1921-1944
    rate2[3,284:307]=354.6
    # 1945-1960
    rate2[3,308:323]=196.3
    # 1961-1972
    rate2[3,324:335]=83.7
    # 1973-1990
    rate2[3,336:353]=39
    # 1991-1997
    rate2[3,354:360]=0
    # Overwriting the default rate
    rate<-rate2[3,1:360]
    
    # Creating a discrete transition rate
    rate_floor <- floor(rate)
    rate_decimals <- (rate-trunc(rate))*10
    
    rounding_random <- runif(1997-(1638+1),0.0,10.0)
    
    for (s in 1:length(rate)){
      rate[s] <- ifelse (rounding_random[s]<rate_decimals[s],(rate_floor[s]+1),rate_floor[s])
    }
  }
  
  # Correction for the rate (the last rate numbers are NA for some undiscovered reason)
  rate[359:360]=0
  
  # Calculating the forest extent according to the rate
  forest_ext <- matrix(ncol = 360,nrow = 1)
  for (i in 1:360){
    if (i==1){
      forest_ext[1]=length(which(forest_1638==2))
    }
    if (i>1){
      forest_ext[i] = forest_ext[i-1]-rate[i] 
    }
  }
  
  ## Dynamic part: the transitions
  ##########################################
  # Loop for time
  for (t in start_year:end_year){
    v=v+1
    
    # Calculating the suitability per period
    if (historical_suitability==1){
      if (t==1638){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1642){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1646){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1647){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1648){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1649){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1650){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1651){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1653){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1655){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1656){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1658){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1664){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1665){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1667){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1668){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1669){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1673){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1675){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_br*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1679){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1683){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1686){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1690){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1693){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1698){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_br*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1703){
        suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
      }
      if (t==1722){
        suitab_total <- (suitab_agr_l*w2_agr_suitab+suitab_dem*w2_dem+suitab_slope*w2_slope+suitab_rivers_1772*w2_riv_french+suitab_sett_1772*w2_sett_french)/(w2_agr_suitab+w2_dem+w2_slope+w2_sett_french+w2_riv_french) 
      } 
      if (t==1810){
        suitab_total <- (suitab_agr_l*w3_agr_suitab+suitab_coast_dist*w3_coast_dist+suitab_dem*w3_dem+suitab_slope*w3_slope)/(w3_agr_suitab+w3_coast_dist+w3_dem+w3_slope)  
      }
      if (t==1835){
        suitab_total <- (suitab_agr_l*w4_agr_suitab+suitab_coast_dist*w4_coast_dist+suitab_dem*w4_dem+suitab_slope*w4_slope)/(w4_agr_suitab+w4_coast_dist+w4_dem+w4_slope)  
      }
      if (t==1872){
        suitab_total <- (suitab_agr_l*w5_agr_suitab+suitab_coast_dist*w5_coast_dist+suitab_dem*w5_dem+suitab_slope*w5_slope)/(w5_agr_suitab+w5_coast_dist+w5_dem+w5_slope)  
      }
      if (t==1968){
        suitab_total <- (suitab_agr_l*w6_agr_suitab+suitab_coast_dist*w6_coast_dist+suitab_dem*w6_dem+suitab_slope*w6_slope)/(w6_agr_suitab+w6_coast_dist+w6_dem+w6_slope)  
      }
    }
    
    # Selecting the forest cells
    forest_indices <- which(forest==2)
    suitab_part <- suitab_total[forest_indices]
    
    # Determine cells to transition
    rate_adj <- length(forest_indices)-forest_ext[v]
    
    transition_number = rate_adj
    suitab_part <- sort(suitab_part,decreasing=TRUE)
    if (transition_number>0){
      suitab_part <- suitab_part[1:transition_number]
    }
    if (transition_number==0){
      suitab_part <-1000000
    }
    
    # Determine cell indices
    cell_indices=0
    
    cell_indices <- which(suitab_total %in% suitab_part)
    
    # Transition cells
    forest[cell_indices]=1
    
    # Print status
    cat("\014")
    print(t)
    forest_record <- append(forest_record,length(which(forest==2)))
    rate_adj_record <- append(rate_adj_record,rate_adj)
    
    # Save calibration years
    if (start_year==1638){
      if (t==1638){
        f_1638 <- rotate(forest)
      }  
    }
    if (start_year<1686 & end_year>1684){
      if (t==1684){
        f_1685 <- rotate(forest)
      }  
    }
    if (start_year<1773 & end_year>1771){
      if (t==1772){
        f_1773 <- rotate(forest)
      }  
    }
    if (start_year<1835 & end_year>1833){
      if (t==1834){
        f_1835 <- rotate(forest)
      }  
    }
    if (start_year<1873 & end_year>1871){
      if (t==1871){
        f_1872 <- rotate(forest)
      }  
    }
    if (start_year<1936 & end_year>1934){
      if (t==1934){
        f_1935 <- rotate(forest)
      }  
    }
    if (start_year<1998 & end_year>1996){
      if (t==1995){
        f_1997 <- rotate(forest)
      }  
    }
    
    # Dynamic visualisation
    if (dynamic_visualisation==1){
      p = t/5-trunc(t/5)
      if (p==0){
        forest1 <- as.matrix(aggregate(ts(forest, frequency=3), 1, sum))
        forest1 <- as.data.frame(t(matrix(forest1, 200)))
        forest1 <- as.matrix(aggregate(ts(forest1, frequency=3), 1, sum))
        forest1 <- t(forest1)
        forest1 <- rotate(forest1)
        graphics.off()
        image(forest1,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
        main_title <- paste("Forest extent in ",toString(t))
        title(main=main_title)
      }  
    }
    
    # Adjusting regrowth (turned off for the Dutch period)
    if (t>1685){
      regrowth[forest==1]=regrowth[forest==1]+1
      forest[regrowth>19]=2
      regrowth_record = append(regrowth_record,length(which(regrowth>19)))
      regrowth[regrowth>19]=0
      regrowth[forest==2]=0  
    }
  }
  
  # Calculate the accuracy of the calibration years
  acc_1685=1-length(which(abs(f_1685-rotate(forest_1685))==1))/cells_island
  acc_1773=1-length(which(abs(f_1773-rotate(forest_1773))==1))/cells_island
  acc_1835=1-length(which(abs(f_1835-rotate(forest_1835))==1))/cells_island
  acc_1872=1-length(which(abs(f_1872-rotate(forest_1872))==1))/cells_island
  acc_1935=1-length(which(abs(f_1935-rotate(forest_1935))==1))/cells_island
  acc_1997=1-length(which(abs(f_1997-rotate(forest_1997))==1))/cells_island
}

if (calibration==1){
  
  for (m in 1:calibration_runs){
    cat("\014")
    print(m)
    
    # Improved rate (more historical)
    if (improved_rate==1){
      rate2 <- matrix(nrow = 3,ncol=360)
      rate2[1,1:360] = 1638:1997
      rate2[2,1:360] = 1:360
      # 1638
      rate2[3,1]=0
      # 1639-1645
      rate2[3,2:8]=175
      # 1646-1650
      rate2[3,9:13]=220
      # 1651-1655
      rate2[3,14:18]=80
      # 1656-1658
      rate2[3,19:21]=118
      # 1659-1663
      rate2[3,22:26]=0
      # 1664-1668
      rate2[3,27:31]=170
      # 1669-1685
      rate2[3,32:47]=145
      # 1686-1705
      rate2[3,48:69]=100
      # 1706-1721
      rate2[3,70:84]=0
      # 1722-1766
      rate2[3,85:129]=479
      # 1766-1790
      rate2[3,130:153]=697.666
      # 1791-1845
      rate2[3,154:208]=905.9
      # 1846-1870
      rate2[3,209:233]=1647.5
      # 1871-1890
      rate2[3,234:253]=1011
      # 1891-1920
      rate2[3,254:283]=281
      # 1921-1944
      rate2[3,284:307]=354.6
      # 1945-1960
      rate2[3,308:323]=196.3
      # 1961-1972
      rate2[3,324:335]=83.7
      # 1973-1990
      rate2[3,336:353]=39
      # 1991-1997
      rate2[3,354:360]=0
      # Overwriting the default rate
      rate<-rate2[3,1:360]
      
      # Creating a discrete transition rate
      rate_floor <- floor(rate)
      rate_decimals <- (rate-trunc(rate))*10
      
      rounding_random <- runif(1997-(1638+1),0.0,10.0)
      
      for (s in 1:length(rate)){
        rate[s] <- ifelse (rounding_random[s]<rate_decimals[s],(rate_floor[s]+1),rate_floor[s])
      }
    }
    
    # Correction for the rate (the last rate numbers are NA for some undiscovered reason)
    rate[359:360]=0
    
    # Calculating the forest extent according to the rate
    forest_ext <- matrix(ncol = 360,nrow = 1)
    for (i in 1:360){
      if (i==1){
        forest_ext[1]=length(which(forest_1638==2))
      }
      if (i>1){
        forest_ext[i] = forest_ext[i-1]-rate[i] 
      }
    }
    
    ## Dynamic part: the transitions
    ##########################################
    # Loop for time
    for (t in start_year:end_year){
      v=v+1
      
      if (t==1638){
        v=1
        forest=forest_1638
      }
      
      # Calculating the suitability per period
      if (historical_suitability==1){
        if (t==1638){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1642){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1646){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1647){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1648){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1649){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1650){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1651){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1653){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1655){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1656){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1658){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1664){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_map*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1665){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1667){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1668){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1669){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1673){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1675){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_br*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1679){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1683){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_bat*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1686){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_gp*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1690){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_ted*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1693){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_flacq*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1698){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_br*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1703){
          suitab_total <- (suitab_agr_l*w1_agr_suitab+suitab_dem*w1_dem+suitab_slope*w1_slope+suitab_pl*w1_sett_dutch)/(w1_agr_suitab+w1_dem+w1_slope+w1_sett_dutch)  
        }
        if (t==1722){
          suitab_total <- (suitab_agr_l*w2_agr_suitab+suitab_dem*w2_dem+suitab_slope*w2_slope+suitab_rivers_1772*w2_riv_french+suitab_sett_1772*w2_sett_french)/(w2_agr_suitab+w2_dem+w2_slope+w2_sett_french+w2_riv_french) 
        } 
        if (t==1810){
          suitab_total <- (suitab_agr_l*w3_agr_suitab+suitab_coast_dist*w3_coast_dist+suitab_dem*w3_dem+suitab_slope*w3_slope)/(w3_agr_suitab+w3_coast_dist+w3_dem+w3_slope)  
        }
        if (t==1835){
          suitab_total <- (suitab_agr_l*w4_agr_suitab+suitab_coast_dist*w4_coast_dist+suitab_dem*w4_dem+suitab_slope*w4_slope)/(w4_agr_suitab+w4_coast_dist+w4_dem+w4_slope)  
        }
        if (t==1872){
          suitab_total <- (suitab_agr_l*w5_agr_suitab+suitab_coast_dist*w5_coast_dist+suitab_dem*w5_dem+suitab_slope*w5_slope)/(w5_agr_suitab+w5_coast_dist+w5_dem+w5_slope)  
        }
        if (t==1968){
          suitab_total <- (suitab_agr_l*w6_agr_suitab+suitab_coast_dist*w6_coast_dist+suitab_dem*w6_dem+suitab_slope*w6_slope)/(w6_agr_suitab+w6_coast_dist+w6_dem+w6_slope)  
        }
      }
      
      # Selecting the forest cells
      forest_indices <- which(forest==2)
      suitab_part <- suitab_total[forest_indices]
      
      # Determine cells to transition
      rate_adj <- length(forest_indices)-forest_ext[v]
      
      transition_number = rate_adj
      suitab_part <- sort(suitab_part,decreasing=TRUE)
      if (transition_number>0){
        suitab_part <- suitab_part[1:transition_number]
      }
      if (transition_number==0){
        suitab_part <-1000000
      }
      
      # Determine cell indices
      cell_indices=0
      
      cell_indices <- which(suitab_total %in% suitab_part)
      
      # Transition cells
      forest[cell_indices]=1
      
      # Print status
      forest_record <- append(forest_record,length(which(forest==2)))
      rate_adj_record <- append(rate_adj_record,rate_adj)
      
      # Save calibration years
      if (start_year==1638){
        if (t==1638){
          f_1638 <- rotate(forest)
        }  
      }
      if (start_year<1686 & end_year>1684){
        if (t==1684){
          f_1685 <- rotate(forest)
        }  
      }
      if (start_year<1773 & end_year>1771){
        if (t==1772){
          f_1773 <- rotate(forest)
        }  
      }
      if (start_year<1835 & end_year>1833){
        if (t==1834){
          f_1835 <- rotate(forest)
        }  
      }
      if (start_year<1873 & end_year>1871){
        if (t==1871){
          f_1872 <- rotate(forest)
        }  
      }
      if (start_year<1936 & end_year>1934){
        if (t==1934){
          f_1935 <- rotate(forest)
        }  
      }
      if (start_year<1998 & end_year>1996){
        if (t==1995){
          f_1997 <- rotate(forest)
        }  
      }
      
      # Adjusting regrowth
      if (t>1685){
        regrowth[forest==1]=regrowth[forest==1]+1
        forest[regrowth>19]=2
        regrowth_record = append(regrowth_record,length(which(regrowth>19)))
        regrowth[regrowth>19]=0
        regrowth[forest==2]=0  
      }
    }
    
    # Calculate the accuracy of the calibration years
    acc_1685=append(acc_1685,1-length(which(abs(f_1685-rotate(forest_1685))==1))/cells_island)
    acc_1773=append(acc_1773,1-length(which(abs(f_1773-rotate(forest_1773))==1))/cells_island)
    acc_1835=append(acc_1835,1-length(which(abs(f_1835-rotate(forest_1835))==1))/cells_island)
    acc_1872=append(acc_1872,1-length(which(abs(f_1872-rotate(forest_1872))==1))/cells_island)
    acc_1935=append(acc_1935,1-length(which(abs(f_1935-rotate(forest_1935))==1))/cells_island)
    acc_1997=append(acc_1997,1-length(which(abs(f_1997-rotate(forest_1997))==1))/cells_island)
    acc_engind=append(acc_engind,(acc_1835[length(acc_1835)]+acc_1872[length(acc_1872)]+acc_1935[length(acc_1935)]+acc_1997[length(acc_1997)])/4)
      
    # Change the weight of the appropriate factor
    if (cal_w1_agr_suitab==1){
      w1_agr_suitab=w1_agr_suitab+calibration_steps  
    }
    if (cal_w1_dem==1){
      w1_dem=w1_dem+calibration_steps  
    }
    if (cal_w1_slope==1){
      w1_slope=w1_slope+calibration_steps  
    }
    if (cal_w1_sett_dutch==1){
      w1_sett_dutch=w1_sett_dutch+calibration_steps  
    }
    if (cal_w2_agr_suitab==1){
      w2_agr_suitab=w2_agr_suitab+calibration_steps
    }
    if (cal_w2_dem==1){
      w2_dem=w2_dem+calibration_steps
    }
    if (cal_w2_slope==1){
      w2_slope=w2_slope+calibration_steps
    }
    if (cal_w2_sett_french==1){
      w2_sett_french=w2_sett_french+calibration_steps
    }
    if (cal_w2_riv_french==1){
      w2_riv_french=w2_riv_french+calibration_steps
    }
    if (cal_w3_agr_suitab==1){
      w3_agr_suitab=w3_agr_suitab+calibration_steps
    }
    if (cal_w3_dem==1){
      w3_dem=w3_dem+calibration_steps
    }
    if (cal_w3_slope==1){
      w3_slope=w3_slope+calibration_steps
    }
    if (cal_w3_coast_dist==1){
      w3_coast_dist=w3_coast_dist+calibration_steps
    }
    if (cal_w4_agr_suitab==1){
      w4_agr_suitab=w4_agr_suitab+calibration_steps
    }
    if (cal_w4_dem==1){
      w4_dem=w4_dem+calibration_steps
    }
    if (cal_w4_slope==1){
      w4_slope=w4_slope+calibration_steps
    }
    if (cal_w4_coast_dist==1){
      w4_coast_dist=w4_coast_dist+calibration_steps
    }
    if (cal_w5_agr_suitab==1){
      w5_agr_suitab=w5_agr_suitab+calibration_steps
    }
    if (cal_w5_dem==1){
      w5_dem=w5_dem+calibration_steps
    }
    if (cal_w5_slope==1){
      w5_slope=w5_slope+calibration_steps
    }
    if (cal_w5_coast_dist==1){
      w5_coast_dist=w5_coast_dist+calibration_steps
    }
    if (cal_w6_agr_suitab==1){
      w6_agr_suitab=w6_agr_suitab+calibration_steps
    }
    if (cal_w6_dem==1){
      w6_dem=w6_dem+calibration_steps
    }
    if (cal_w6_slope==1){
      w6_slope=w6_slope+calibration_steps
    }
    if (cal_w6_coast_dist==1){
      w6_coast_dist=w6_coast_dist+calibration_steps
    }
  }
}

# Correcting arrays
forest_record=forest_record[2:length(forest_record)]
rate_adj_record=rate_adj_record[2:length(rate_adj_record)]
regrowth_record=regrowth_record[2:length(regrowth_record)]

## Visualisations
#########################################################

# Deforestation rate
if (graph_deforestation_rate==1){
  plot(start_year:end_year,rate,type="l",xlab = "Years",ylab = "Deforestation rate (hectares/year)")
  title("Deforestation rate on Mauritius")
}
  
# Forest extent
if (graph_forest_extent==1){
  plot(start_year:end_year,forest_record,type="l",xlab = "Years",ylab = "Forest extent (hectares)")
  title("Forest extent on Mauritius")
}

# Images of forest extent
if (images_forest_extent==1){  
  par(mfrow=c(2,3))
  image(x,y,f_1685,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1685")
  image(x,y,f_1773,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1773")
  image(x,y,f_1835,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1835")
  image(x,y,f_1872,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1872")
  image(x,y,f_1935,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1935")
  image(x,y,f_1997,axes=FALSE,xaxt="n",yaxt="n",col=mycols,xlab = "",ylab = "")
  title(main="Forest extent in 1997")
}

# Regrowth
if (graph_regrowth==1){
  plot(1686:end_year,regrowth_record,type="l",xlab = "Years",ylab = "Regrowth (hectares)")
  title("Regrowth (secondary forest) on Mauritius")
}

# Create graph of the calibration array
if (cal_w1_agr_suitab==1){
  plot(acc_1685[2:length(acc_1685)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the Dutch occupation")
}
if (cal_w1_dem==1){
  plot(acc_1685[2:length(acc_1685)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the Dutch occupation")
}
if (cal_w1_slope==1){
  plot(acc_1685[2:length(acc_1685)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during the Dutch occupation")
}
if (cal_w1_sett_dutch==1){
  plot(acc_1685[2:length(acc_1685)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the settlement factor during the Dutch occupation")
}

if (cal_w2_agr_suitab==1){
  plot(acc_1773[2:length(acc_1773)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the French occupation")
}
if (cal_w2_dem==1){
  plot(acc_1773[2:length(acc_1773)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the French occupation")
}
if (cal_w2_slope==1){
  plot(acc_1773[2:length(acc_1773)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during the French occupation")
}
if (cal_w2_sett_french==1){
  plot(acc_1773[2:length(acc_1773)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the settlements factor during the French occupation")
}
if (cal_w2_riv_french==1){
  plot(acc_1773[2:length(acc_1773)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the rivers factor during the French occupation")
}

if (cal_w3_agr_suitab==1){
  plot(acc_1835[2:length(acc_1835)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the early period of the English occupation")
}
if (cal_w3_dem==1){
  plot(acc_1835[2:length(acc_1835)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the early period of the English occupation")
}
if (cal_w3_slope==1){
  plot(acc_1835[2:length(acc_1835)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during the early period of the English occupation")
}
if (cal_w3_coast_dist==1){
  plot(acc_1835[2:length(acc_1835)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the coast distance factor during the early period of the English occupation")
}

if (cal_w4_agr_suitab==1){
  plot(acc_1872[2:length(acc_1872)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the middle period of the English occupation")
}
if (cal_w4_dem==1){
  plot(acc_1872[2:length(acc_1872)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the middle period of the English occupation")
}
if (cal_w4_slope==1){
  plot(acc_1872[2:length(acc_1872)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during middle period of the English occupation")
}
if (cal_w4_coast_dist==1){
  plot(acc_1872[2:length(acc_1872)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the coast distance factor during the middle period of the English occupation")
}

if (cal_w5_agr_suitab==1){
  plot(acc_1935[2:length(acc_1935)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the late period of the English occupation")
}
if (cal_w5_dem==1){
  plot(acc_1935[2:length(acc_1935)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the late period of the English occupation")
}
if (cal_w5_slope==1){
  plot(acc_1935[2:length(acc_1935)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during the late period of the English occupation")
}
if (cal_w5_coast_dist==1){
  plot(acc_1935[2:length(acc_1935)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the coast distance factor during the late period of the English occupation and independence")
}

if (cal_w6_agr_suitab==1){
  plot(acc_1997[2:length(acc_1997)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the agricultural suitability factor during the independence")
}
if (cal_w6_dem==1){
  plot(acc_1997[2:length(acc_1997)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the DEM factor during the independence")
}
if (cal_w6_slope==1){
  plot(acc_1997[2:length(acc_1997)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the slope factor during the independence")
}
if (cal_w6_coast_dist==1){
  plot(acc_1997[2:length(acc_1997)],type="l",xlab = "Weight",ylab = "Accuracy (fraction)")
  title("Accuracy of the coast distance factor during the independence")
}

# Check correct rate
if (check_rate==1){
  print(length(which(forest_1685==2))-length(which(f_1685==2)))
  print(length(which(forest_1773==2))-length(which(f_1773==2)))
  print(length(which(forest_1835==2))-length(which(f_1835==2)))
  print(length(which(forest_1872==2))-length(which(f_1872==2)))
  print(length(which(forest_1935==2))-length(which(f_1935==2)))
  print(length(which(forest_1997==2))-length(which(f_1997==2)))
}

# Print final accuracy
if (print_accuracy==1){
  print(signif(acc_1685[length(acc_1685)]),digits=4)
  print(signif(acc_1773[length(acc_1773)]),digits=4)
  print(signif(acc_1835[length(acc_1835)]),digits=4)
  print(signif(acc_1872[length(acc_1872)]),digits=4)
  print(signif(acc_1935[length(acc_1935)]),digits=4)
  print(signif(acc_1997[length(acc_1997)]),digits=4)
  print(signif(acc_1835[length(acc_1835)]+acc_1872[length(acc_1872)]+acc_1935[length(acc_1935)]+acc_1997[length(acc_1997)])/4)
}