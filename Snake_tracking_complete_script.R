###############################################################################
###############################################################################
###############################################################################
###########################                         ###########################
###########################                         ###########################
###########################  Snake tracking script  ###########################
###########################                         ###########################
###########################                         ###########################
###############################################################################
###############################################################################
###############################################################################

# Part 0: Packages
##### 

library(viridis)
library(readxl)
library(writexl)
library(ggplot2)
library(tidyr)
library(plyr)
library(sp)
library(lubridate)
library(dplyr)
library(irr)
library(utils)
library(tibble)
library(EnvStats)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(Hmisc)
library(ggpubr)
library(vegan)
library(sm)
library(survminer)
library(ggfortify)
library(factoextra)
library(graphics)

#####

# Part 1: Load required data
#####

# Set working directory
setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking")

# Load file with snake info
Snake_data <- read_excel("Snake_captures.xlsx")

# Subset tracked snakes
Snake_data_tracked <- Snake_data[Snake_data$Tracking %in% c("Y"), ]

# Subset tracked + genetics/isotopes + x-ray sample 
Snake_data_tracked <- Snake_data_tracked[!is.na(Snake_data_tracked$SVL), ]


# Load snake 3D array
load("Snake_array.Rda")

# Load arena 3D array
load("Arena_array.Rda")

# Load summary file from all snake trials
summary_tracking <- read_excel("summary_tracking.xlsx")

# Load summary file from all snake trials + environmental variables
summary_tracking_complete <- read_excel("summary_tracking_complete.xlsx")

#####

# Part 2: Data extraction for each individual snake and 3D storing (already done)
#####

# Variable names for each column
variable_names <- c("t", "x", "y", "Central_area", "Intermediate_area", "Refuge_area_1", "Refuge_area_2", "Initial_refuge", "Distant_wall", "Lateral_1", "Lateral_2", "Arena_area", "Area", "dist_hole_0", "dist_hole_1", "dist_hole_2", "dist_hole_3", "dist_hole_4", "dist_min", "hole_dist_min", "dist_traveled", "behaviour")

# We create empty data frame that is going to store: x axis = time step, y = variables, z = individuals
Snake_array <-  array(NA, dim = c(3600, length(variable_names), nrow(Snake_data_tracked)))

# We set the column names for the variables across the array
dimnames(Snake_array)[[2]] <- variable_names

# We set the ID code for all tracked snakes across the array
dimnames(Snake_array)[[3]] <- Snake_data_tracked$Snake_ID


# For each tracked individual 
# THIS TAKES A LOT OF TIME. DO NOT RUN AND LOAD DIRECTLY FROM WD
for(i in 1:nrow(Snake_data_tracked)){
  
  # If it is the first trial (no repeatability)...
  if(Snake_data_tracked$Repeatability[i] == "N"){
    
    # ... set the working directory in the correct folder...
    setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking/Excels tracking")
    
    # ... and if it didn't get out of the original refuge, fill the correspondent data layer with 0. Otherwise...
    if(Snake_data_tracked$Time_head_out_sec_raw[i] == "NA"){
      
      # We copy an already complete file
      coordinates_complete <- as.data.frame(Snake_array[, , 1])
      
      # And fill with 0
      coordinates_complete$x <- 0
      coordinates_complete$y <- 0
      coordinates_complete$Central_area <- 0
      coordinates_complete$Intermediate_area <- 0
      coordinates_complete$Refuge_area_1 <- 0
      coordinates_complete$Refuge_area_2 <- 0
      coordinates_complete$Initial_refuge <- 0
      coordinates_complete$Distant_wall <- 0
      coordinates_complete$Lateral_1 <- 0
      coordinates_complete$Lateral_2 <- 0
      coordinates_complete$Arena_area <- 0 
      coordinates_complete$Area <- "Initial refuge"
      coordinates_complete$dist_hole_0 <- 0
      coordinates_complete$dist_hole_1 <- 0.86
      coordinates_complete$dist_hole_2 <- 1.3
      coordinates_complete$dist_hole_3 <- 1.3
      coordinates_complete$dist_hole_4 <- 0.86
      coordinates_complete$dist_min <- 0
      coordinates_complete$hole_dist_min <- "dist_hole_0"
      coordinates_complete$dist_traveled <- 0
      coordinates_complete$behaviour <- "hidden"
      
      # It might happen that the snake is being recorded for more than 1h, and we want to know their behaviour in 1h-tests. It's for this reason that we are going to cut all videos to 1h (3600 sec). If the video lasts less than 1h, for the moment we are not going to do anything.
      if(nrow(coordinates_complete) > 3600){
        
        coordinates_complete <- coordinates_complete[1:3600, ]
        
      } else {
        
      }
      
      Snake_array[1:nrow(coordinates_complete), , i] <- as.matrix(coordinates_complete[ , ])
      
      assign("Snake_array", Snake_array, envir=globalenv())
      
      paste("Snake_array", i)
      
      
      # ... extract data from its tracked Excel file.
    } else if(Snake_data_tracked$Time_head_out_sec_raw[i] != "NA"){
      
      # Select data based on snake's ID
      data_name <- paste0(Snake_data_tracked$Snake_ID[i], "_coordinates.xlsx")
      
      # Standardized name for data
      coordinates <- read_excel(data_name)
      
      # Select arena limits and refuges from specific individual
      name_arena <- paste0(Snake_data_tracked$Snake_ID[i], "_arena.xlsx")
      
      arena <- read_excel(name_arena)
      
      # How much time does the snake hide before first head-out
      sec_bf_head_out <- as.numeric(Snake_data_tracked$Time_head_out_sec_raw[i])
      
      # Let's find which seconds is the snake missing in the recording (hidden)
      missing_time <- min(coordinates$t):max(coordinates$t)
      
      # How many seconds does the snake stay inside the refuges after coming out the initial refuge
      missing_time_recorded <- length(missing_time[!missing_time %in% coordinates$t])
      
      # Which seconds is the snake missing
      missing_time_sec <- missing_time[!missing_time %in% coordinates$t]
      
      
      # See in which area is the snake at each time step
      
      # The snake, even if it is outside the arena, is still over a certain area, so we are going to widen the limits of the areas for detecting in which area is the snake even if it is outside the arena limits. It does not have an effect over the behaviour
      
      off_limit <- 0.2
      
      # Central area
      Central_area <- point.in.polygon(coordinates$x, 
                                       coordinates$y, 
                                       c(arena[1, 2]- 0.75, arena[2, 2] - 0.75, arena[3, 2] + 0.75, arena[4, 2] + 0.75), 
                                       c(arena[1, 3] + 0.5, arena[2, 3] - 0.5, arena[3, 3] - 0.5, arena[4, 3] + 0.5))
      
      coordinates$Central_area <- Central_area
      
      # Intermediate area
      Intermediate_area <- point.in.polygon(coordinates$x, 
                                            coordinates$y, 
                                            c(arena[1, 2]- 0.5, arena[2, 2] - 0.5, arena[3, 2] + 0.5, arena[4, 2] + 0.5), 
                                            c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25))
      
      coordinates$Intermediate_area <- Intermediate_area
      
      coordinates$Intermediate_area <- coordinates$Intermediate_area - coordinates$Central_area
      
      # Refuge area 1
      Refuge_area_1 <- point.in.polygon(coordinates$x, 
                                        coordinates$y, 
                                        c(arena[1, 2]- 0.25, arena[2, 2] - 0.25, arena[2, 2] - 0.5, arena[2, 2] - 0.5),
                                        c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[2, 3] - 0.25, arena[1, 3] + 0.25))
      
      coordinates$Refuge_area_1 <- Refuge_area_1
      
      # Refuge area 2
      Refuge_area_2 <- point.in.polygon(coordinates$x, 
                                        coordinates$y, 
                                        c(arena[4, 2] + 0.5, arena[3, 2] + 0.5, arena[3, 2] + 0.25, arena[4, 2] + 0.25),
                                        c(arena[4, 3] + 0.25, arena[3, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25))
      
      coordinates$Refuge_area_2 <- Refuge_area_2
      
      
      # Initial refuge area
      Initial_refuge <- point.in.polygon(coordinates$x, 
                                         coordinates$y, 
                                         c(arena[1, 2] - 0.25, arena[1, 2] -0.25, arena[4, 2] + 0.25, arena[4, 2] + 0.25),
                                         c(arena[1, 3] - off_limit, arena[1, 3] + 0.25, arena[4, 3] + 0.25, arena[4, 3] - off_limit))
      
      coordinates$Initial_refuge <- Initial_refuge
      
      
      # Distant wall area
      Distant_wall <- point.in.polygon(coordinates$x, 
                                       coordinates$y, 
                                       c(arena[2, 2] - 0.25, arena[2, 2] - 0.25, arena[3, 2] + 0.25, arena[3, 2] + 0.25),
                                       c(arena[2, 3] - 0.25, arena[2, 3] + off_limit, arena[3, 3] + off_limit, arena[3, 3] - 0.25))
      
      coordinates$Distant_wall <- Distant_wall
      
      
      # Lateral refuge area 1
      Lateral_1 <- point.in.polygon(coordinates$x, 
                                    coordinates$y, 
                                    c(arena[1, 2] + off_limit, arena[2, 2] + off_limit, arena[2, 2] - 0.25, arena[1, 2] - 0.25),
                                    c(arena[1, 3] - off_limit, arena[2, 3] + off_limit, arena[2, 3] + off_limit, arena[1, 3] - off_limit))
      
      coordinates$Lateral_1 <- Lateral_1
      
      # Lateral refuge area 2
      Lateral_2 <- point.in.polygon(coordinates$x, 
                                    coordinates$y, 
                                    c(arena[3, 2] - off_limit, arena[4, 2] - off_limit, arena[4, 2] + 0.25, arena[3, 2] + 0.25),
                                    c(arena[3, 3] + off_limit, arena[4, 3] - off_limit, arena[4, 3] - off_limit, arena[3, 3] + off_limit))
      
      coordinates$Lateral_2 <- Lateral_2
      
      
      # Arena area. We are going to expand the arena limits 2cm on each side, in order to take into account edge exploration as inside the arena, not ouside
      Arena_area <- point.in.polygon(coordinates$x, 
                                     coordinates$y, 
                                     c(arena[1, 2] + 0.02, arena[2, 2] + 0.02, arena[3, 2] - 0.02, arena[4, 2] - 0.02),
                                     c(arena[1, 3] - 0.02, arena[2, 3] + 0.02, arena[3, 3] + 0.02, arena[4, 3] - 0.02))
      
      coordinates$Arena_area <- Arena_area
      
      # Column that indicates the area the snake is at every time step
      for(j in 1:nrow(coordinates)){
        
        if(coordinates$Arena_area[j] == 1 | coordinates$Arena_area[j] == 0){
          
          if(coordinates$Lateral_2[j] == 1){
            
            coordinates$Area[j] <- "Lateral 2"
            
          } else if(coordinates$Lateral_1[j] == 1){
            
            coordinates$Area[j] <- "Lateral 1"
            
          } else if(coordinates$Distant_wall[j] == 1){
            
            coordinates$Area[j] <- "Distant wall"
            
          } else if(coordinates$Refuge_area_1[j] == 1){
            
            coordinates$Area[j] <- "Refuge area 1"
            
          } else if(coordinates$Refuge_area_2[j] == 1) {
            
            coordinates$Area[j] <- "Refuge area 2"
            
          } else if(coordinates$Initial_refuge[j] == 1) {
            
            coordinates$Area[j] <- "Initial refuge"
            
          } else if(coordinates$Intermediate_area[j] == 1) {
            
            coordinates$Area[j] <- "Intermediate area"
            
          } else if(coordinates$Central_area[j] == 1){
            
            coordinates$Area[j] <- "Central area"
            
          } else  {
            
            if (j == 1){
              
              # Initial frame
              coordinates$Area[j] <- "Initial refuge"
              
            } else {
              
              # If any of the conditions are met, repeat previous area
              coordinates$Area[j] <- coordinates$Area[j - 1]
              
            }
            
          }
          
        }
        
      }
      
      
      
      # Distance to each hole
      
      # Entrance hole (hole_0)
      hole_0 <- c(arena[5, 2], arena[5, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_0[k] <- sqrt((coordinates$x[k] - as.numeric(hole_0[1]))^2 + (coordinates$y[k] - as.numeric(hole_0[2]))^2)
        
      }
      
      # Hole 1
      hole_1 <- c(arena[6, 2], arena[6, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_1[k] <- sqrt((coordinates$x[k] - as.numeric(hole_1[1]))^2 + (coordinates$y[k] - as.numeric(hole_1[2]))^2)
        
      }
      
      # Hole 2
      hole_2 <- c(arena[7, 2], arena[7, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_2[k] <- sqrt((coordinates$x[k] - as.numeric(hole_2[1]))^2 + (coordinates$y[k] - as.numeric(hole_2[2]))^2)
        
      }
      
      # Hole 3
      hole_3 <- c(arena[8, 2], arena[8, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_3[k] <- sqrt((coordinates$x[k] - as.numeric(hole_3[1]))^2 + (coordinates$y[k] - as.numeric(hole_3[2]))^2)
        
      }
      
      # Hole 4
      hole_4 <- c(arena[9, 2], arena[9, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_4[k] <- sqrt((coordinates$x[k] - as.numeric(hole_4[1]))^2 + (coordinates$y[k] - as.numeric(hole_4[2]))^2)
        
      }
      
      
      # Which is the hole the snake is closest to at each moment
      for(k in 1:nrow(coordinates)){
        
        a <- coordinates$dist_hole_0[k]
        b <- coordinates$dist_hole_1[k]
        c <- coordinates$dist_hole_2[k]
        d <- coordinates$dist_hole_3[k]
        e <- coordinates$dist_hole_4[k]
        
        # Write the distance to the closest hole
        coordinates$dist_min[k] <- min(c(a, b, c, d, e))
        
        # And the name of that hole
        coordinates$hole_dist_min[k] <- colnames(coordinates[which((coordinates[k, ] == coordinates$dist_min[k]) == TRUE)])[1]
        
      }
      
      
      # Which distance has the snake traveled from 1 sec to the next one
      dist_traveled <- data.frame(matrix(ncol = 1, nrow = nrow(coordinates)))
      
      for(k in 2:nrow(coordinates)){
        
        dist_traveled[k, 1] <- sqrt((coordinates$x[k] - coordinates$x[k - 1])^2 + (coordinates$y[k] - coordinates$y[k - 1])^2)
        
      }
      
      dist_traveled[1, 1] <- 0
      
      coordinates$dist_traveled <- as.numeric(dist_traveled[, 1])
      
      
      # Which behaviour is the snake exhibiting?
      for(k in 1:nrow(coordinates)){
        
        if(coordinates$Arena_area[k] == 0){
          
          if(coordinates$dist_traveled[k] < 0.01){
            
            coordinates$behaviour[k] <- "hiding / head-out"
            
          } else {
            
            coordinates$behaviour[k] <- "exploring wall"
            
          }
          
        } else {
          
          if(coordinates$dist_traveled[k] < 0.01 && coordinates$dist_min[k] < 0.1){
            
            coordinates$behaviour[k] <- "head-out"
            
          } else if(coordinates$dist_traveled[k] < 0.01 && coordinates$dist_min[k] > 0.1){
            
            coordinates$behaviour[k] <- "immobile exposed"
            
          } else {
            
            coordinates$behaviour[k] <- "moving"
            
          }
          
        }
        
      }
      
      
      # We put the correct time steps to each frame
      coordinates$t <- coordinates$t + as.numeric(sec_bf_head_out) + 1
      
      # We create a data set that is going to keep all the information to what happens during the whole video
      
      # How long is the total video (nº sec recorded + nº sec snake is hidden bf exiting for the first time + nº seconds hidden during recording)
      total_frames <- as.numeric(sec_bf_head_out) + nrow(coordinates) + missing_time_recorded
      
      # Column names
      col_names_coordinates <- colnames(coordinates)
      
      # We create the data frame
      coordinates_complete <- data.frame(matrix(ncol = ncol(coordinates), nrow = total_frames))
      
      # And put the column names
      colnames(coordinates_complete) <- col_names_coordinates
      
      # Seconds that the whole recordind lasts
      coordinates_complete$t <- seq(1:total_frames)
      
      # Now, we assign each recorded parameter to its correspondent second on the 1h record
      for(k in 1:nrow(coordinates)){
        
        values_coordinates <- coordinates[k, ]
        
        coordinates_complete[as.numeric(values_coordinates[1]), ] <- values_coordinates
        
        
      }
      
      
      # Values for when the snake is still inside the initial refuge
      for(j in 1:nrow(coordinates_complete)){
        
        # If the snake hasn't exited the initial refuge
        if(j <= as.numeric(sec_bf_head_out)){
          
          coordinates_complete$x[j] <- 0
          coordinates_complete$y[j] <- 0
          coordinates_complete$Central_area[j] <- 0
          coordinates_complete$Intermediate_area[j] <- 0
          coordinates_complete$Refuge_area_1[j] <- 0
          coordinates_complete$Refuge_area_2[j] <- 0
          coordinates_complete$Initial_refuge[j] <- 0
          coordinates_complete$Distant_wall[j] <- 0
          coordinates_complete$Lateral_1[j] <- 0
          coordinates_complete$Lateral_2[j] <- 0
          coordinates_complete$Arena_area[j] <- 0 
          coordinates_complete$Area[j] <- "Initial refuge"
          coordinates_complete$dist_hole_0[j] <- 0
          coordinates_complete$dist_hole_1[j] <- coordinates_complete$dist_hole_1[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_2[j] <- coordinates_complete$dist_hole_2[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_3[j] <- coordinates_complete$dist_hole_3[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_4[j] <- coordinates_complete$dist_hole_4[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_min[j] <- 0
          coordinates_complete$hole_dist_min[j] <- "dist_hole_0"
          coordinates_complete$dist_traveled[j] <- 0
          coordinates_complete$behaviour[j] <- "hidden"
          
          
        } else {
          
          # If the snake is missing (hidden), which parameters are going to exhibit those frames. We are going to copy the parameters from the last recorded frame, except specific parameters.
          if(is.na(coordinates_complete$x[j]) == TRUE) {
            
            coordinates_complete$x[j] <- coordinates_complete$x[j - 1]
            coordinates_complete$y[j] <- coordinates_complete$y[j - 1]
            coordinates_complete$Central_area[j] <- coordinates_complete$Central_area[j - 1]
            coordinates_complete$Intermediate_area[j] <- coordinates_complete$Intermediate_area[j - 1]
            coordinates_complete$Refuge_area_1[j] <- coordinates_complete$Refuge_area_1[j - 1]
            coordinates_complete$Refuge_area_2[j] <- coordinates_complete$Refuge_area_2[j - 1]
            coordinates_complete$Initial_refuge[j] <-  coordinates_complete$Initial_refuge[j - 1]
            coordinates_complete$Distant_wall[j] <- coordinates_complete$Distant_wall[j - 1]
            coordinates_complete$Lateral_1[j] <- coordinates_complete$Lateral_1[j - 1]
            coordinates_complete$Lateral_2[j] <- coordinates_complete$Lateral_2[j - 1]
            coordinates_complete$Arena_area[j] <- coordinates_complete$Arena_area[j - 1]
            coordinates_complete$Area[j] <- coordinates_complete$Area[j - 1]
            coordinates_complete$dist_hole_0[j] <- coordinates_complete$dist_hole_0[j - 1]
            coordinates_complete$dist_hole_1[j] <- coordinates_complete$dist_hole_1[j - 1]
            coordinates_complete$dist_hole_2[j] <- coordinates_complete$dist_hole_2[j - 1]
            coordinates_complete$dist_hole_3[j] <- coordinates_complete$dist_hole_3[j - 1]
            coordinates_complete$dist_hole_4[j] <- coordinates_complete$dist_hole_4[j - 1]
            coordinates_complete$dist_min[j] <- 0
            coordinates_complete$hole_dist_min[j] <- coordinates_complete$hole_dist_min[j - 1]
            coordinates_complete$dist_traveled[j] <- 0
            coordinates_complete$behaviour[j] <- "hidden"
            
          }
          
        }
        
      }
      
      # It might happen that the snake is being recorded for more than 1h, and we want to know their behaviour in 1h-tests. It's for this reason that we are going to cut all videos to 1h (3600 sec). If the video lasts less than 1h, for the moment we are not going to do anything.
      if(nrow(coordinates_complete) > 3600){
        
        coordinates_complete <- coordinates_complete[1:3600, ]
        
      } else {
        
      }
      
      Snake_array[1:nrow(coordinates_complete), , i] <- as.matrix(coordinates_complete[ , ])
      
      assign("Snake_array", Snake_array, envir=globalenv())
      
      paste("Snake_array", i)
      
    }
    
    # If it is a snake used for repeatability, ...
  } else if(Snake_data_tracked$Repeatability[i] == "Y"){
    
    # ... set the correct wd ...
    setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking/Excels tracking/Replicates")
    
    # ... and check if the snake did get out of the initial refuge or not. If it didn't, fill with 0.
    if(Snake_data_tracked$Time_head_out_sec_raw[i] == "NA"){
      
      
      # We copy an already complete file
      coordinates_complete <- as.data.frame(Snake_array[, , 1])
      
      # And fill with 0
      coordinates_complete$x <- 0
      coordinates_complete$y <- 0
      coordinates_complete$Central_area <- 0
      coordinates_complete$Intermediate_area <- 0
      coordinates_complete$Refuge_area_1 <- 0
      coordinates_complete$Refuge_area_2 <- 0
      coordinates_complete$Initial_refuge <- 0
      coordinates_complete$Distant_wall <- 0
      coordinates_complete$Lateral_1 <- 0
      coordinates_complete$Lateral_2 <- 0
      coordinates_complete$Arena_area <- 0 
      coordinates_complete$Area <- "Initial refuge"
      coordinates_complete$dist_hole_0 <- 0
      coordinates_complete$dist_hole_1 <- 0.86
      coordinates_complete$dist_hole_2 <- 1.3
      coordinates_complete$dist_hole_3 <- 1.3
      coordinates_complete$dist_hole_4 <- 0.86
      coordinates_complete$dist_min <- 0
      coordinates_complete$hole_dist_min <- "dist_hole_0"
      coordinates_complete$dist_traveled <- 0
      coordinates_complete$behaviour <- "hidden"
      
      # It might happen that the snake is being recorded for more than 1h, and we want to know their behaviour in 1h-tests. It's for this reason that we are going to cut all videos to 1h (3600 sec). If the video lasts less than 1h, for the moment we are not going to do anything.
      if(nrow(coordinates_complete) > 3600){
        
        coordinates_complete <- coordinates_complete[1:3600, ]
        
      } else {
        
      }
      
      Snake_array[1:nrow(coordinates_complete), , i] <- as.matrix(coordinates_complete[ , ])
      
      assign("Snake_array", Snake_array, envir=globalenv())
      
      paste("Snake_array", i)
      
      
      
      # If the snake did get out of the refuge, extract data from excel
    } else if(Snake_data_tracked$Time_head_out_sec_raw[i] != "NA"){
      
      # Select data based on snake's ID
      data_name <- paste0(Snake_data_tracked$Snake_ID[i], "_coordinates.xlsx")
      
      # Standardized name for data
      coordinates <- read_excel(data_name)
      
      # Select arena limits and refuges from specific individual
      name_arena <- paste0(Snake_data_tracked$Snake_ID[i], "_arena.xlsx")
      
      arena <- read_excel(name_arena)
      
      # How much time does the snake hide before first head-out
      sec_bf_head_out <- as.numeric(Snake_data_tracked$Time_head_out_sec_raw[i])
      
      # Let's find which seconds is the snake missing in the recording (hidden)
      missing_time <- min(coordinates$t):max(coordinates$t)
      
      # How many seconds does the snake stay inside the refuges after coming out the initial refuge
      missing_time_recorded <- length(missing_time[!missing_time %in% coordinates$t])
      
      # Which seconds is the snake missing
      missing_time_sec <- missing_time[!missing_time %in% coordinates$t]
      
      
      # See in which area is the snake at each time step
      
      # The snake, even if it is outside the arena, is still over a certain area, so we are going to widen the limits of the areas for detecting in which area is the snake even if it is outside the arena limits. It does not have an effect over the behaviour
      
      off_limit <- 0.2
      
      # Central area
      Central_area <- point.in.polygon(coordinates$x, 
                                       coordinates$y, 
                                       c(arena[1, 2]- 0.75, arena[2, 2] - 0.75, arena[3, 2] + 0.75, arena[4, 2] + 0.75), 
                                       c(arena[1, 3] + 0.5, arena[2, 3] - 0.5, arena[3, 3] - 0.5, arena[4, 3] + 0.5))
      
      coordinates$Central_area <- Central_area
      
      # Intermediate area
      Intermediate_area <- point.in.polygon(coordinates$x, 
                                            coordinates$y, 
                                            c(arena[1, 2]- 0.5, arena[2, 2] - 0.5, arena[3, 2] + 0.5, arena[4, 2] + 0.5), 
                                            c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25))
      
      coordinates$Intermediate_area <- Intermediate_area
      
      coordinates$Intermediate_area <- coordinates$Intermediate_area - coordinates$Central_area
      
      # Refuge area 1
      Refuge_area_1 <- point.in.polygon(coordinates$x, 
                                        coordinates$y, 
                                        c(arena[1, 2]- 0.25, arena[2, 2] - 0.25, arena[2, 2] - 0.5, arena[2, 2] - 0.5),
                                        c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[2, 3] - 0.25, arena[1, 3] + 0.25))
      
      coordinates$Refuge_area_1 <- Refuge_area_1
      
      # Refuge area 2
      Refuge_area_2 <- point.in.polygon(coordinates$x, 
                                        coordinates$y, 
                                        c(arena[4, 2] + 0.5, arena[3, 2] + 0.5, arena[3, 2] + 0.25, arena[4, 2] + 0.25),
                                        c(arena[4, 3] + 0.25, arena[3, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25))
      
      coordinates$Refuge_area_2 <- Refuge_area_2
      
      
      # Initial refuge area
      Initial_refuge <- point.in.polygon(coordinates$x, 
                                         coordinates$y, 
                                         c(arena[1, 2] - 0.25, arena[1, 2] -0.25, arena[4, 2] + 0.25, arena[4, 2] + 0.25),
                                         c(arena[1, 3] - off_limit, arena[1, 3] + 0.25, arena[4, 3] + 0.25, arena[4, 3] - off_limit))
      
      coordinates$Initial_refuge <- Initial_refuge
      
      
      # Distant wall area
      Distant_wall <- point.in.polygon(coordinates$x, 
                                       coordinates$y, 
                                       c(arena[2, 2] - 0.25, arena[2, 2] - 0.25, arena[3, 2] + 0.25, arena[3, 2] + 0.25),
                                       c(arena[2, 3] - 0.25, arena[2, 3] + off_limit, arena[3, 3] + off_limit, arena[3, 3] - 0.25))
      
      coordinates$Distant_wall <- Distant_wall
      
      
      # Lateral refuge area 1
      Lateral_1 <- point.in.polygon(coordinates$x, 
                                    coordinates$y, 
                                    c(arena[1, 2] + off_limit, arena[2, 2] + off_limit, arena[2, 2] - 0.25, arena[1, 2] - 0.25),
                                    c(arena[1, 3] - off_limit, arena[2, 3] + off_limit, arena[2, 3] + off_limit, arena[1, 3] - off_limit))
      
      coordinates$Lateral_1 <- Lateral_1
      
      # Lateral refuge area 2
      Lateral_2 <- point.in.polygon(coordinates$x, 
                                    coordinates$y, 
                                    c(arena[3, 2] - off_limit, arena[4, 2] - off_limit, arena[4, 2] + 0.25, arena[3, 2] + 0.25),
                                    c(arena[3, 3] + off_limit, arena[4, 3] - off_limit, arena[4, 3] - off_limit, arena[3, 3] + off_limit))
      
      coordinates$Lateral_2 <- Lateral_2
      
      
      # Arena area. We are going to expand the arena limits 2cm on each side, in order to take into account edge exploration as inside the arena, not ouside
      Arena_area <- point.in.polygon(coordinates$x, 
                                     coordinates$y, 
                                     c(arena[1, 2] + 0.02, arena[2, 2] + 0.02, arena[3, 2] - 0.02, arena[4, 2] - 0.02),
                                     c(arena[1, 3] - 0.02, arena[2, 3] + 0.02, arena[3, 3] + 0.02, arena[4, 3] - 0.02))
      
      coordinates$Arena_area <- Arena_area
      
      # Column that indicates the area the snake is at every time step
      for(j in 1:nrow(coordinates)){
        
        if(coordinates$Arena_area[j] == 1 | coordinates$Arena_area[j] == 0){
          
          if(coordinates$Lateral_2[j] == 1){
            
            coordinates$Area[j] <- "Lateral 2"
            
          } else if(coordinates$Lateral_1[j] == 1){
            
            coordinates$Area[j] <- "Lateral 1"
            
          } else if(coordinates$Distant_wall[j] == 1){
            
            coordinates$Area[j] <- "Distant wall"
            
          } else if(coordinates$Refuge_area_1[j] == 1){
            
            coordinates$Area[j] <- "Refuge area 1"
            
          } else if(coordinates$Refuge_area_2[j] == 1) {
            
            coordinates$Area[j] <- "Refuge area 2"
            
          } else if(coordinates$Initial_refuge[j] == 1) {
            
            coordinates$Area[j] <- "Initial refuge"
            
          } else if(coordinates$Intermediate_area[j] == 1) {
            
            coordinates$Area[j] <- "Intermediate area"
            
          } else if(coordinates$Central_area[j] == 1){
            
            coordinates$Area[j] <- "Central area"
            
          } else  {
            
            if (j == 1){
              
              # Initial frame
              coordinates$Area[j] <- "Initial refuge"
              
            } else {
              
              # If any of the conditions are met, repeat previous area
              coordinates$Area[j] <- coordinates$Area[j - 1]
              
            }
            
          }
          
        }
        
      }
      
      
      
      # Distance to each hole
      
      # Entrance hole (hole_0)
      hole_0 <- c(arena[5, 2], arena[5, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_0[k] <- sqrt((coordinates$x[k] - as.numeric(hole_0[1]))^2 + (coordinates$y[k] - as.numeric(hole_0[2]))^2)
        
      }
      
      # Hole 1
      hole_1 <- c(arena[6, 2], arena[6, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_1[k] <- sqrt((coordinates$x[k] - as.numeric(hole_1[1]))^2 + (coordinates$y[k] - as.numeric(hole_1[2]))^2)
        
      }
      
      # Hole 2
      hole_2 <- c(arena[7, 2], arena[7, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_2[k] <- sqrt((coordinates$x[k] - as.numeric(hole_2[1]))^2 + (coordinates$y[k] - as.numeric(hole_2[2]))^2)
        
      }
      
      # Hole 3
      hole_3 <- c(arena[8, 2], arena[8, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_3[k] <- sqrt((coordinates$x[k] - as.numeric(hole_3[1]))^2 + (coordinates$y[k] - as.numeric(hole_3[2]))^2)
        
      }
      
      # Hole 4
      hole_4 <- c(arena[9, 2], arena[9, 3])
      
      for(k in 1:nrow(coordinates)){
        
        coordinates$dist_hole_4[k] <- sqrt((coordinates$x[k] - as.numeric(hole_4[1]))^2 + (coordinates$y[k] - as.numeric(hole_4[2]))^2)
        
      }
      
      
      # Which is the hole the snake is closest to at each moment
      for(k in 1:nrow(coordinates)){
        
        a <- coordinates$dist_hole_0[k]
        b <- coordinates$dist_hole_1[k]
        c <- coordinates$dist_hole_2[k]
        d <- coordinates$dist_hole_3[k]
        e <- coordinates$dist_hole_4[k]
        
        # Write the distance to the closest hole
        coordinates$dist_min[k] <- min(c(a, b, c, d, e))
        
        # And the name of that hole
        coordinates$hole_dist_min[k] <- colnames(coordinates[which((coordinates[k, ] == coordinates$dist_min[k]) == TRUE)])[1]
        
      }
      
      
      # Which distance has the snake traveled from 1 sec to the next one
      dist_traveled <- data.frame(matrix(ncol = 1, nrow = nrow(coordinates)))
      
      for(k in 2:nrow(coordinates)){
        
        dist_traveled[k, 1] <- sqrt((coordinates$x[k] - coordinates$x[k - 1])^2 + (coordinates$y[k] - coordinates$y[k - 1])^2)
        
      }
      
      dist_traveled[1, 1] <- 0
      
      coordinates$dist_traveled <- as.numeric(dist_traveled[, 1])
      
      
      # Which behaviour is the snake exhibiting?
      for(k in 1:nrow(coordinates)){
        
        if(coordinates$Arena_area[k] == 0){
          
          if(coordinates$dist_traveled[k] < 0.01){
            
            coordinates$behaviour[k] <- "hiding / head-out"
            
          } else {
            
            coordinates$behaviour[k] <- "exploring wall"
            
          }
          
        } else {
          
          if(coordinates$dist_traveled[k] < 0.01 && coordinates$dist_min[k] < 0.1){
            
            coordinates$behaviour[k] <- "head-out"
            
          } else if(coordinates$dist_traveled[k] < 0.01 && coordinates$dist_min[k] > 0.1){
            
            coordinates$behaviour[k] <- "immobile exposed"
            
          } else {
            
            coordinates$behaviour[k] <- "moving"
            
          }
          
        }
        
      }
      
      
      # We put the correct time steps to each frame
      coordinates$t <- coordinates$t + as.numeric(sec_bf_head_out) + 1
      
      # We create a data set that is going to keep all the information to what happens during the whole video
      
      # How long is the total video (nº sec recorded + nº sec snake is hidden bf exiting for the first time + nº seconds hidden during recording)
      total_frames <- as.numeric(sec_bf_head_out) + nrow(coordinates) + missing_time_recorded
      
      # Column names
      col_names_coordinates <- colnames(coordinates)
      
      # We create the data frame
      coordinates_complete <- data.frame(matrix(ncol = ncol(coordinates), nrow = total_frames))
      
      # And put the column names
      colnames(coordinates_complete) <- col_names_coordinates
      
      # Seconds that the whole recordind lasts
      coordinates_complete$t <- seq(1:total_frames)
      
      # Now, we assign each recorded parameter to its correspondent second on the 1h record
      for(k in 1:nrow(coordinates)){
        
        values_coordinates <- coordinates[k, ]
        
        coordinates_complete[as.numeric(values_coordinates[1]), ] <- values_coordinates
        
        
      }
      
      
      # Values for when the snake is still inside the initial refuge
      for(j in 1:nrow(coordinates_complete)){
        
        # If the snake hasn't exited the initial refuge
        if(j <= as.numeric(sec_bf_head_out)){
          
          coordinates_complete$x[j] <- 0
          coordinates_complete$y[j] <- 0
          coordinates_complete$Central_area[j] <- 0
          coordinates_complete$Intermediate_area[j] <- 0
          coordinates_complete$Refuge_area_1[j] <- 0
          coordinates_complete$Refuge_area_2[j] <- 0
          coordinates_complete$Initial_refuge[j] <- 0
          coordinates_complete$Distant_wall[j] <- 0
          coordinates_complete$Lateral_1[j] <- 0
          coordinates_complete$Lateral_2[j] <- 0
          coordinates_complete$Arena_area[j] <- 0 
          coordinates_complete$Area[j] <- "Initial refuge"
          coordinates_complete$dist_hole_0[j] <- 0
          coordinates_complete$dist_hole_1[j] <- coordinates_complete$dist_hole_1[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_2[j] <- coordinates_complete$dist_hole_2[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_3[j] <- coordinates_complete$dist_hole_3[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_hole_4[j] <- coordinates_complete$dist_hole_4[as.numeric(sec_bf_head_out) + 1]
          coordinates_complete$dist_min[j] <- 0
          coordinates_complete$hole_dist_min[j] <- "dist_hole_0"
          coordinates_complete$dist_traveled[j] <- 0
          coordinates_complete$behaviour[j] <- "hidden"
          
          
        } else {
          
          # If the snake is missing (hidden), which parameters are going to exhibit those frames. We are going to copy the parameters from the last recorded frame, except specific parameters.
          if(is.na(coordinates_complete$x[j]) == TRUE) {
            
            coordinates_complete$x[j] <- coordinates_complete$x[j - 1]
            coordinates_complete$y[j] <- coordinates_complete$y[j - 1]
            coordinates_complete$Central_area[j] <- coordinates_complete$Central_area[j - 1]
            coordinates_complete$Intermediate_area[j] <- coordinates_complete$Intermediate_area[j - 1]
            coordinates_complete$Refuge_area_1[j] <- coordinates_complete$Refuge_area_1[j - 1]
            coordinates_complete$Refuge_area_2[j] <- coordinates_complete$Refuge_area_2[j - 1]
            coordinates_complete$Initial_refuge[j] <-  coordinates_complete$Initial_refuge[j - 1]
            coordinates_complete$Distant_wall[j] <- coordinates_complete$Distant_wall[j - 1]
            coordinates_complete$Lateral_1[j] <- coordinates_complete$Lateral_1[j - 1]
            coordinates_complete$Lateral_2[j] <- coordinates_complete$Lateral_2[j - 1]
            coordinates_complete$Arena_area[j] <- coordinates_complete$Arena_area[j - 1]
            coordinates_complete$Area[j] <- coordinates_complete$Area[j - 1]
            coordinates_complete$dist_hole_0[j] <- coordinates_complete$dist_hole_0[j - 1]
            coordinates_complete$dist_hole_1[j] <- coordinates_complete$dist_hole_1[j - 1]
            coordinates_complete$dist_hole_2[j] <- coordinates_complete$dist_hole_2[j - 1]
            coordinates_complete$dist_hole_3[j] <- coordinates_complete$dist_hole_3[j - 1]
            coordinates_complete$dist_hole_4[j] <- coordinates_complete$dist_hole_4[j - 1]
            coordinates_complete$dist_min[j] <- 0
            coordinates_complete$hole_dist_min[j] <- coordinates_complete$hole_dist_min[j - 1]
            coordinates_complete$dist_traveled[j] <- 0
            coordinates_complete$behaviour[j] <- "hidden"
            
          }
          
        }
        
      }
      
      # It might happen that the snake is being recorded for more than 1h, and we want to know their behaviour in 1h-tests. It's for this reason that we are going to cut all videos to 1h (3600 sec). If the video lasts less than 1h, for the moment we are not going to do anything.
      if(nrow(coordinates_complete) > 3600){
        
        coordinates_complete <- coordinates_complete[1:3600, ]
        
      } else {
        
      }
      
      Snake_array[1:nrow(coordinates_complete), , i] <- as.matrix(coordinates_complete[ , ])
      
      assign("Snake_array", Snake_array, envir=globalenv())
      
      paste("Snake_array", i)
      
    }
    
  }
  
}

# Set WD to save 3D array
setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking")

# Save 3D array
save(Snake_array, file="Snake_array.Rda")

# Load snake 3D array
load("Snake_array.Rda")


# Arena array
# Variable names for each column
variable_names <- c("t", "x", "y")

# We create empty data frame that is going to store: x axis = time step, y = variables, z = individuals
Arena_array <- array(NA, dim = c(13, length(variable_names), nrow(Snake_data_tracked)))

# We set the column names for the variables across the array
dimnames(Arena_array)[[2]] <- variable_names

# We set the ID code for all tracked snakes across the array
dimnames(Arena_array)[[3]] <- Snake_data_tracked$Snake_ID

# Arena array
# THIS TAKES A LOT OF TIME. DO NOT RUN AND LOAD DIRECTLY FROM WD
for(i in 1:nrow(Snake_data_tracked)){

  # If it is the first trial (no repeatability)...
  if(Snake_data_tracked$Repeatability[i] == "N"){
    
    # ... set the working directory in the correct folder...
    setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking/Excels tracking")
    
    # ... and if it didn't get out of the original refuge, fill the correspondent data layer with 0. Otherwise...
    if(Snake_data_tracked$Time_head_out_sec_raw[i] == "NA"){
      
      # There is no arena drown for snakes that didn't get out of the refuge, so we have to draw them manually
      arena <- as.data.frame(matrix(NA, nrow = 13, ncol = 3))
      
      colnames(arena) <- c("t", "x", "y")
      
      arena$t <- c(seq(0, 12, 1))
      
      arena$x <- c(1, 1, -1, -1, 0, 0.73, 0.73, -0.73, -0.73, 0.73, 0.73, -0.73, -0.73)
      
      arena$y <- c(0, 1.5, 1.5, 0, 0, 0.43, 1.06, 1.06, 0.43, 0, 1.5, 1.5, 0)
      
      Arena_array[1:nrow(arena), , i] <- as.matrix(arena[ , ])
      
      assign("Arena_array", Arena_array, envir=globalenv())
      
      paste("Arena_array", i)
      
      
      # ... extract data from its tracked Excel file.
    } else if(Snake_data_tracked$Time_head_out_sec_raw[i] != "NA"){
      
      # Select arena limits and refuges from specific individual
      name_arena <- paste0(Snake_data_tracked$Snake_ID[i], "_arena.xlsx")
      
      arena <- read_excel(name_arena)
      
      Arena_array[1:nrow(arena), , i] <- as.matrix(arena[ , ])
      
      assign("Arena_array", Arena_array, envir=globalenv())
      
      paste("Arena_array", i)
      
    }
    
    # If it is a snake used for repeatability, ...
  } else if(Snake_data_tracked$Repeatability[i] == "Y"){
    
    # ... set the correct wd ...
    setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking/Excels tracking/Replicates")
    
    # ... and check if the snake did get out of the initial refuge or not. If it didn't, fill with 0.
    if(Snake_data_tracked$Time_head_out_sec_raw[i] == "NA"){
      
      # There is no arena drown for snakes that didn't get out of the refuge, so we have to draw them manually
      arena <- as.data.frame(matrix(NA, nrow = 13, ncol = 3))
      
      colnames(arena) <- c("t", "x", "y")
      
      arena$t <- c(seq(0, 12, 1))
      
      arena$x <- c(1, 1, -1, -1, 0, 0.73, 0.73, -0.73, -0.73, 0.73, 0.73, -0.73, -0.73)
      
      arena$y <- c(0, 1.5, 1.5, 0, 0, 0.43, 1.06, 1.06, 0.43, 0, 1.5, 1.5, 0)
      
      Arena_array[1:nrow(arena), , i] <- as.matrix(arena[ , ])
      
      assign("Arena_array", Arena_array, envir=globalenv())
      
      paste("Arena_array", i)
      
      
      # If the snake did get out of the refuge, extract data from excel
    } else if(Snake_data_tracked$Time_head_out_sec_raw[i] != "NA"){
      
      # Select arena limits and refuges from specific individual
      name_arena <- paste0(Snake_data_tracked$Snake_ID[i], "_arena.xlsx")
      
      arena <- read_excel(name_arena)
      
      Arena_array[1:nrow(arena), , i] <- as.matrix(arena[ , ])
      
      assign("Arena_array", Arena_array, envir=globalenv())
      
      paste("Arena_array", i)
      
    }
    
  }
  
  
  
}

# Set WD to save 3D array
setwd("C:/Users/marc9/Desktop/Marc/CREAF/Snake tracking")

# Save 3D array
save(Arena_array, file="Arena_array.Rda")

# Load arena 3D array
load("Arena_array.Rda")


#####


# Part 3: Individual plots
#####

# ID: Snake ID
# repeatability: If it was the first trial of the snake ("N") or it was a repeatability trial ("Y")
# time_limit: how many seconds do we want to analyse from each snake

ID <- 14739
repeatability <- "N"
time_limit <- 3000

data_visualization <- function(ID, repeatability, time_limit){
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  # Extract arena info from array
  arena <- as.data.frame(Arena_array[, , selected_ind])
  
  # Let's start drawing
  plot(c(-1.1, 1.1), c(-0.1, 1.6), col = "white", xlab = "X", ylab = "Y") 
  
  # Lateral refuge area 1
  polygon(x = c(arena[1, 2], arena[2, 2], arena[2, 2] - 0.25, arena[1, 2] - 0.25),
          y = c(arena[1, 3], arena[2, 3], arena[2, 3], arena[1, 3]),
          col = viridis(6)[6], border = viridis(6)[6]) 
  
  # Lateral refuge area 2
  polygon(x = c(arena[3, 2], arena[4, 2], arena[4, 2] + 0.25, arena[3, 2] + 0.25),
          y = c(arena[3, 3], arena[4, 3], arena[4, 3], arena[3, 3]),
          col = viridis(6)[6], border = viridis(6)[6])
  
  # Distant wall area
  polygon(x = c(arena[2, 2] - 0.25, arena[2, 2] - 0.25, arena[3, 2] + 0.25, arena[3, 2] + 0.25),
          y = c(arena[2, 3] - 0.25, arena[2, 3], arena[3, 3], arena[3, 3] - 0.25),
          col = viridis(6)[5], border = viridis(6)[5]) 
  
  # Initial refuge area
  polygon(x =c(arena[1, 2] - 0.25, arena[1, 2] -0.25, arena[4, 2] + 0.25, arena[4, 2] + 0.25),
          y = c(arena[1, 3], arena[1, 3] + 0.25, arena[4, 3] + 0.25, arena[4, 3]),
          col = viridis(6)[4], border = viridis(6)[4]) 
  
  # Refuge area 1
  polygon(x = c(arena[1, 2]- 0.25, arena[2, 2] - 0.25, arena[2, 2] - 0.5, arena[1, 2] - 0.5),
          y = c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[2, 3] - 0.25, arena[1, 3] + 0.25),
          col = viridis(6)[3], border = viridis(6)[3]) 
  
  # Refuge area 2
  polygon(x = c(arena[4, 2] + 0.5, arena[3, 2] + 0.5, arena[3, 2] + 0.25, arena[4, 2] + 0.25),
          y = c(arena[4, 3] + 0.25, arena[3, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25),
          col = viridis(6)[3], border = viridis(6)[3])
  
  # Intermediate area
  polygon(x = c(arena[1, 2]- 0.5, arena[2, 2] - 0.5, arena[3, 2] + 0.5, arena[4, 2] + 0.5),
          y = c(arena[1, 3] + 0.25, arena[2, 3] - 0.25, arena[3, 3] - 0.25, arena[4, 3] + 0.25),
          col = viridis(6)[2], border = viridis(6)[2]) 
  
  # Central area
  polygon(x = c(arena[1, 2]- 0.75, arena[2, 2] - 0.75, arena[3, 2] + 0.75, arena[4, 2] + 0.75),
          y = c(arena[1, 3] + 0.5, arena[2, 3] - 0.5, arena[3, 3] - 0.5, arena[4, 3] + 0.5),
          col = viridis(6)[1], border = viridis(6)[1]) 
  
  
  # Border real arena
  polygon(x = c(arena[1, 2], arena[2, 2], arena[3, 2], arena[4, 2]),
          y = c(arena[1, 3], arena[2, 3], arena[3, 3], arena[4, 3]),  
          lwd = 2) 
  
  # Refuge coordinates
  refuges <- data.frame(cbind(c(arena[5, 2], arena[6, 2], arena[7, 2], arena[8, 2], arena[9, 2]),
                              c(arena[5, 3], arena[6, 3], arena[7, 3], arena[8, 3], arena[9, 3])))
  
  # Plot refuges
  points(refuges, pch = 19, cex = 4, col = "white")
  
  # Plot snakee trail
  points(coordinates$x, coordinates$y, type ="l", col = "black", lwd = 2)
  
  # In case you want to draw lateral refuge limits 
  # lines(c(arena[10, 2], arena[11, 2]), c(arena[10, 3], arena[11, 3]))
  
  # lines(c(arena[12, 2], arena[13, 2]), c(arena[12, 3], arena[13, 3]))
  
  # Exit from refuge (Last frame you were missing, and now I have detected you)
  for(i in 2:nrow(coordinates)){
    
    if(coordinates$behaviour[i - 1] == "hidden" && coordinates$behaviour[i] != "hidden"){
      
      points(coordinates$x[i], coordinates$y[i], 
             cex = 1, pch = 24, bg = "white", col = "black", lwd = 2)
      
    } else {
      
      
    }
    
  }
  
  
  # Enter to refuge (Now I see you, but next frame I don't)
  for(i in 1:(nrow(coordinates) - 1)){
    
    if(coordinates$behaviour[i] != "hidden" && coordinates$behaviour[i + 1] == "hidden"){
      
      points(coordinates$x[i], coordinates$y[i], 
             cex = 1, pch = 25, bg = "white", col = "black", lwd = 2)
      
    } else {
      
      
    }
    
  }
  
}

density_function <- function(ID, repeatability, time_limit){
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  no_hidden <- coordinates[!grepl("hidden", coordinates$behaviour), ]
  
  coordinates <- no_hidden
  
  ggplot(coordinates, aes(x=as.numeric(x), y=as.numeric(y))) +
    geom_bin2d(bins = 20) +
    scale_fill_continuous(type = "viridis") +
    labs(x = "x", y = "y") +
    theme_bw()
  
  # ggplot(coordinates, aes(x=x, y=y) ) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")
  
}

graph_area <- function(ID, repeatability, time_limit){
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  # How many seconds the snake has the snake been in each region
  
  Time_central_area <- 0
  
  Time_intermediate_area <- 0
  
  Time_refuge_area_1 <- 0
  
  Time_refuge_area_2 <- 0
  
  Time_initial_refuge <- 0
  
  Time_distant_wall <- 0
  
  Time_lateral_area_1 <- 0
  
  Time_lateral_area_2 <- 0
  
  Time_hidden <- 0
  
  # How much time does the snake spend in each area
  for(i in 1:nrow(coordinates)){
    
    if(coordinates$behaviour[i] == "hidden"){
      
      Time_hidden <- Time_hidden + 1
      
    } else {
      
      if (coordinates$Area[i] == "Central area") {
        
        Time_central_area <- Time_central_area + 1
        
      } else if (coordinates$Area[i] == "Intermediate area") {
        
        Time_intermediate_area <- Time_intermediate_area + 1
        
      } else if (coordinates$Area[i] == "Refuge area 1"){
        
        Time_refuge_area_1 <- Time_refuge_area_1 + 1
        
      } else if (coordinates$Area[i] == "Refuge area 2"){
        
        Time_refuge_area_2 <- Time_refuge_area_2 + 1
        
      } else if (coordinates$Area[i] == "Initial refuge"){
        
        Time_initial_refuge <- Time_initial_refuge + 1
        
      } else if (coordinates$Area[i] == "Distant wall"){
        
        Time_distant_wall <- Time_distant_wall + 1
        
      } else if (coordinates$Area[i] == "Lateral 1"){
        
        Time_lateral_area_1 <- Time_lateral_area_1 + 1
        
      } else if (coordinates$Area[i] == "Lateral 2"){
        
        Time_lateral_area_2 <- Time_lateral_area_2 + 1
        
      }
      
    }
    
  }
  
  
  # Vector with number of seconds
  total_use  <- c(Time_hidden, Time_lateral_area_1, Time_lateral_area_2, Time_distant_wall, Time_initial_refuge, Time_refuge_area_1, Time_refuge_area_2, Time_intermediate_area, Time_central_area)
  # Vector with percentage 
  percent_total_use <- round(total_use/sum(total_use)*100, 1)
  
  percent_total_use <- as.data.frame(percent_total_use)
  
  area_names <- c("Hidden", "Lateral refuge 1", "Lateral refuge 2", "Distant wall", "Initial refuge",  "Refuge 1", "Refuge 2", "Intermediate",  "Central")
  
  ggplot(data = percent_total_use, aes(x = area_names, y = percent_total_use)) +
    geom_bar(stat="identity", fill = c("grey50", viridis(6)[6], viridis(6)[6], viridis(6)[5], viridis(6)[4], viridis(6)[3], viridis(6)[3], viridis(6)[2], viridis(6)[1])) + 
    geom_text(aes(label = percent_total_use), vjust = -0.3, size = 3.5) +
    scale_x_discrete(limits = area_names) + 
    labs(title= paste0("Percentage of time spent by snake ", ID , " in each area (", time_limit, " sec)"), 
         x="Area", y = "Percentage (%)") + 
    theme_classic()
  
}

graph_hole <- function(ID, repeatability, time_limit){
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  # How much time does the snake hide before first head-out
  sec_bf_head_out <- as.numeric(Snake_data_tracked$Time_head_out_sec_raw[selected_ind])
  
  # We are not interested in the proximity to any hole when the snake is inside a hole, hidden, so we delete the records when the snake is hidden
  no_hidden <- coordinates[!grepl("hidden", coordinates$behaviour), ]
  
  coordinates <- no_hidden
  
  # Which is the closest hole
  hole_proximity <- table(coordinates$hole_dist_min)
  
  percent_hole_proximity <- round(hole_proximity/sum(hole_proximity)*100, 1)
  
  percent_hole_proximity <- as.data.frame(percent_hole_proximity)
  
  percent_hole_proximity$Var1 <- c("Entrance", "Hole 1", "Hole 2", "Hole 3", "Hole 4")
  
  ggplot(data = percent_hole_proximity, aes(x = Var1, y = Freq)) +
    geom_bar(stat="identity", fill = c(viridis(5)[5], viridis(5)[4], viridis(5)[3], viridis(5)[2], viridis(5)[1])) + 
    geom_text(aes(label = Freq), vjust = -0.3, size = 3.5) +
    scale_x_discrete(limits = percent_hole_proximity$Var1) + 
    labs(title=paste0("Which is the closest refuge while the snake ", ID, " is out? (", time_limit, " sec)"), 
         x="Hole", y = "Percentage (%)") + 
    theme_classic()
  
}

graph_behaviour <- function(ID, repeatability, time_limit){
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  # Which is the behaviour
  moving <- 0
  
  exploring_wall <- 0
  
  head_out <- 0
  
  hidden_head_out <- 0
  
  immobile <- 0
  
  hidden <- 0
  
  for(i in 1:nrow(coordinates)){
    
    if(coordinates$behaviour[i] == "moving"){
      
      moving <- moving + 1
      
    } else if(coordinates$behaviour[i] == "exploring wall"){
      
      exploring_wall <- exploring_wall + 1
      
    } else if(coordinates$behaviour[i] == "head-out"){
      
      head_out <- head_out + 1
      
    } else if(coordinates$behaviour[i] == "hiding / head-out"){
      
      hidden_head_out <- hidden_head_out + 1
      
    } else if(coordinates$behaviour[i] == "immobile exposed"){
      
      immobile <- immobile + 1
      
    } else if(coordinates$behaviour[i] == "hidden"){
      
      hidden <- hidden + 1
      
    }
    
  }
  
  behaviour  <- c(hidden, moving, exploring_wall, head_out, hidden_head_out, immobile)
  
  percent_behaviour <- round(behaviour/sum(behaviour)*100, 1)
  
  percent_behaviour <- as.data.frame(percent_behaviour)
  
  behaviour_names <- c("Hidden", "Moving", "Exploring wall", "Head-out", "Hidden / Head-out", "Immobile exposed")
  
  ggplot(data = percent_behaviour, aes(x = behaviour_names, y = percent_behaviour)) +
    geom_bar(stat="identity", fill = c("grey50", viridis(5)[5], viridis(5)[4], viridis(5)[3], viridis(5)[2], viridis(5)[1])) + 
    geom_text(aes(label = percent_behaviour), vjust = -0.3, size = 3.5) +
    scale_x_discrete(limits = behaviour_names) + 
    labs(title=paste0("Which behaviour is the snake ", ID, " exhibiting? (", time_limit, " sec)"), 
         x="Behaviour", y = "Percentage (%)") + 
    theme_classic()
  
}

speed_plot <- function(ID, repeatability, time_limit) {
  
  par(mfrow = c(2, 1))
  
  par(mar = c(4.5, 5, 0.8, 0.5))  
  
  # Which snake are we going to visualize
  selected_ind <- which(Snake_data_tracked$Snake_ID == ID & Snake_data_tracked$Repeatability == repeatability)
  
  # Extract snake info from array
  coordinates <- as.data.frame(Snake_array[1:time_limit, , selected_ind])
  
  # We are going to delete those recordings where the snake comes out from a refuge, because the speed that the snake might exhibit is not real, it comes from going from refuge A to B from "one frame to the next one"
  for(i in 2:nrow(coordinates)){
    
    if(coordinates$behaviour[i - 1] == "hidden" && coordinates$behaviour[i] != "hidden"){
      
      coordinates$dist_traveled[i] <- 0
      
    }
    
  }
  
  # We plot the results
  plot(c(0, nrow(coordinates)), c(-0.1, max(as.numeric(coordinates$dist_traveled))), 
       col = "white", 
       xlab = "Time (s)", 
       ylab = "Speed (m/s)",
       xaxt = "n", yaxt = "n",
       axes = FALSE)
  
  ticks_x <- seq(0, nrow(coordinates), by = nrow(coordinates)/6)
  
  axis(1, at = ticks_x, 
       labels = c(ticks_x[1],
                  round_any(ticks_x[2], 5),
                  round_any(ticks_x[3], 5),
                  round_any(ticks_x[4], 5),
                  round_any(ticks_x[5], 5),
                  round_any(ticks_x[6], 5),
                  round_any(ticks_x[7], 5)))
  
  ticks_y <- round(max(as.numeric(coordinates$dist_traveled)), 2)
  
  
  axis(2, at = c(0, round((ticks_y/3), 2), round(((ticks_y/3)*2), 2), ticks_y), las = 1)
  
  
  
  # We plot the behaviour underneath
  for(i in 1:nrow(coordinates)){
    
    if(coordinates$behaviour[i] == "hidden") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = "grey50", border = "grey50")
      
    } else if(coordinates$behaviour[i] == "moving") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[5], border = viridis(5)[5])
      
    } else if(coordinates$behaviour[i] == "exploring wall") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[4], border = viridis(5)[4])
      
    } else if(coordinates$behaviour[i] == "head-out") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[3], border = viridis(5)[3])
      
    } else if(coordinates$behaviour[i] == "hiding / head-out") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[2], border = viridis(5)[2])
      
    } else if(coordinates$behaviour[i] == "immobile exposed") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[1], border = viridis(5)[1])
      
    }
    
  }
  
  # We plot the speed
  points(coordinates$t, coordinates$dist_traveled, type = "l")
  
  
  # Now for only the part where the snake has exited the initial refuge
  
  # How much time does the snake hide before first head-out
  sec_bf_head_out <- as.numeric(Snake_data_tracked$Time_head_out_sec_raw[selected_ind])
  
  coordinates_2 <- coordinates[(sec_bf_head_out + 1):nrow(coordinates), ]
  
  plot(c(0, nrow(coordinates_2)), c(-0.1, max(as.numeric(coordinates_2$dist_traveled))),
       col = "white", 
       xlab = "Time (s)", 
       ylab = "Speed (m/s)",
       yaxt = "n", xaxt = "n",
       axes = FALSE
  )
  
  ticks_x <- seq(0, nrow(coordinates_2), by = nrow(coordinates_2)/6)
  
  axis(1, at = ticks_x, 
       labels = c(ticks_x[1] + sec_bf_head_out + 1,
                  round_any(ticks_x[2] + sec_bf_head_out, 5),
                  round_any(ticks_x[3] + sec_bf_head_out, 5),
                  round_any(ticks_x[4] + sec_bf_head_out, 5),
                  round_any(ticks_x[5] + sec_bf_head_out, 5),
                  round_any(ticks_x[6] + sec_bf_head_out, 5),
                  round_any(ticks_x[7] + sec_bf_head_out, 5)))
  
  
  
  
  ticks_y <- round(max(as.numeric(coordinates_2$dist_traveled)), 2)
  
  axis(2, at = c(0, round((ticks_y/3), 2), round(((ticks_y/3)*2), 2), ticks_y), las = 1)
  
  # We plot the behaviour underneath
  for(i in 1:nrow(coordinates_2)){
    
    if(coordinates_2$behaviour[i] == "hidden") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = "grey50", border = "grey50")
      
    } else if(coordinates_2$behaviour[i] == "moving") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[5], border = viridis(5)[5])
      
    } else if(coordinates_2$behaviour[i] == "exploring wall") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[4], border = viridis(5)[4])
      
    } else if(coordinates_2$behaviour[i] == "head-out") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[3], border = viridis(5)[3])
      
    } else if(coordinates_2$behaviour[i] == "hiding / head-out") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[2], border = viridis(5)[2])
      
    } else if(coordinates_2$behaviour[i] == "immobile exposed") {
      
      polygon(x = c((i - 1), (i), (i), (i - 1)),
              y = c(- 0.02, -0.02, -0.08, -0.08),
              col = viridis(5)[1], border = viridis(5)[1])
      
    }
    
  }
  
  # We plot the speed
  points(1:nrow(coordinates_2), coordinates_2$dist_traveled, type = "l")
  
  par(mfrow = c(1, 1))
  
  
}

data_visualization(ID, repeatability, time_limit)

density_function(ID, repeatability, time_limit)

graph_area(ID, repeatability, time_limit)

graph_hole(ID, repeatability, time_limit)

graph_behaviour(ID, repeatability, time_limit)

speed_plot(ID, repeatability, time_limit)

#####


# Part 4: Summarized data storing (already done)
#####

# To create a summary table for all snakes
# This table is going to be maintained for the rest of the analysis

variables_track <- c("ID", "time_central_area", "time_intermediate_area", "time_refuge_area_1", "time_refuge_area_2", "time_initial_refuge", "time_distant_wall", "time_lateral_1", "time_lateral_2", "time_hidden_bf_head_out", "total_time_hidden", "time_entrance", "time_hole_1", "time_hole_2", "time_hole_3", "time_hole_4", "mean_dist_refuge", "mean_dist_closest_refuge", "beh_moving", "beh_exploring_wall", "beh_immobile_exposed", "beh_head-out", "beh_hiding-head-out", "beh_head-out_after_body_out", "avg_speed", "avg_speed_exploring", "total_dist_travelled", "dist_travelled_exploring", "dist_trav_5min_body_out", "dist_trav_10min_body_out", "dist_trav_15min_body_out", "dist_trav_20min_body_out", "dist_trav_25min_body_out", "dist_trav_30min_body_out")

# We create the data frame that is going to store all the data
summary_tracking <- data.frame(matrix(ncol = length(variables_track), nrow = nrow(Snake_data_tracked)))

colnames(summary_tracking) <- variables_track

# We generate the final data frame based on a fixed number of seconds
generate_summary_tracking <- function(time_limit){
  
  for(i in 1:nrow(Snake_data_tracked)){
    
    if(Snake_data_tracked$Time_head_out_sec_raw[i] != "NA" & as.numeric(Snake_data_tracked$Time_head_out_sec_raw[i]) < time_limit){
      
      # Extract snake info from array
      coordinates <- as.data.frame(Snake_array[1:time_limit, , i])
      
      # Extract arena info from array
      arena <- as.data.frame(Arena_array[, , i])
      
      # How much time does the snake hide before first head-out
      sec_bf_head_out <- as.numeric(Snake_data_tracked$Time_head_out_sec_raw[i])
      
      # for those individuals that do head-out, but not body-out
      if(Snake_data_tracked$Time_body_out_sec_raw[i] != "NA"){
        
        sec_bf_body_out <- as.numeric(Snake_data_tracked$Time_body_out_sec_raw[i])
        
      } else  {
        
        sec_bf_body_out <- time_limit
        
      }
      
      
      
      # Let's compile all the info we have
      
      tracking_ind <- summary_tracking[1, ]
      
      tracking_ind[1, ] <- NA
      
      # ID
      tracking_ind$ID <- Snake_data_tracked$Snake_ID[i]
      
      # Area use
      Time_central_area <- 0
      
      Time_intermediate_area <- 0
      
      Time_refuge_area_1 <- 0
      
      Time_refuge_area_2 <- 0
      
      Time_initial_refuge <- 0
      
      Time_distant_wall <- 0
      
      Time_lateral_area_1 <- 0
      
      Time_lateral_area_2 <- 0
      
      Time_hidden <- 0
      
      # How much time does the snake spend in each area
      for(k in 1:nrow(coordinates)){
        
        if(coordinates$behaviour[k] == "hidden"){
          
          Time_hidden <- Time_hidden + 1
          
        } else {
          
          if (coordinates$Area[k] == "Central area") {
            
            Time_central_area <- Time_central_area + 1
            
          } else if (coordinates$Area[k] == "Intermediate area") {
            
            Time_intermediate_area <- Time_intermediate_area + 1
            
          } else if (coordinates$Area[k] == "Refuge area 1"){
            
            Time_refuge_area_1 <- Time_refuge_area_1 + 1
            
          } else if (coordinates$Area[k] == "Refuge area 2"){
            
            Time_refuge_area_2 <- Time_refuge_area_2 + 1
            
          } else if (coordinates$Area[k] == "Initial refuge"){
            
            Time_initial_refuge <- Time_initial_refuge + 1
            
          } else if (coordinates$Area[k] == "Distant wall"){
            
            Time_distant_wall <- Time_distant_wall + 1
            
          } else if (coordinates$Area[k] == "Lateral 1"){
            
            Time_lateral_area_1 <- Time_lateral_area_1 + 1
            
          } else if (coordinates$Area[k] == "Lateral 2"){
            
            Time_lateral_area_2 <- Time_lateral_area_2 + 1
            
          }
          
        }
        
      }
      
      tracking_ind$time_central_area <- Time_central_area
      
      tracking_ind$time_intermediate_area <- Time_intermediate_area
      
      tracking_ind$time_refuge_area_1 <- Time_refuge_area_1
      
      tracking_ind$time_refuge_area_2 <- Time_refuge_area_2
      
      tracking_ind$time_initial_refuge <- Time_initial_refuge
      
      tracking_ind$time_distant_wall <- Time_distant_wall
      
      tracking_ind$time_lateral_1 <- Time_lateral_area_1
      
      tracking_ind$time_lateral_2 <- Time_lateral_area_2
      
      tracking_ind$time_hidden_bf_head_out <- as.numeric(sec_bf_head_out)
      
      tracking_ind$total_time_hidden <- Time_hidden
      
      # Proximity to refuges
      no_hidden <- coordinates[!grepl("hidden", coordinates$behaviour), ]
      
      hole_proximity <- table(no_hidden$hole_dist_min)
      
      tracking_ind$time_entrance <- hole_proximity[1]
      
      tracking_ind$time_hole_1 <- hole_proximity[2]
      
      tracking_ind$time_hole_2 <- hole_proximity[3]
      
      tracking_ind$time_hole_3 <- hole_proximity[4]
      
      tracking_ind$time_hole_4 <- hole_proximity[5]
      
      dist_refuges <- c(no_hidden$dist_hole_0, no_hidden$dist_hole_1, no_hidden$dist_hole_2, no_hidden$dist_hole_3, no_hidden$dist_hole_4)
      
      tracking_ind$mean_dist_refuge <- mean(as.numeric(dist_refuges))
      
      tracking_ind$mean_dist_closest_refuge <- mean(as.numeric(no_hidden$dist_min))
      
      # Behaviour
      moving <- 0
      
      exploring_wall <- 0
      
      head_out <- 0
      
      hidden_head_out <- 0
      
      immobile <- 0
      
      hidden <- 0
      
      
      for(k in 1:nrow(coordinates)){
        
        if(coordinates$behaviour[k] == "moving"){
          
          moving <- moving + 1
          
        } else if(coordinates$behaviour[k] == "exploring wall"){
          
          exploring_wall <- exploring_wall + 1
          
        } else if(coordinates$behaviour[k] == "head-out"){
          
          head_out <- head_out + 1
          
        } else if(coordinates$behaviour[k] == "hiding / head-out"){
          
          hidden_head_out <- hidden_head_out + 1
          
        } else if(coordinates$behaviour[k] == "immobile exposed"){
          
          immobile <- immobile + 1
          
        } else if(coordinates$behaviour[k] == "hidden"){
          
          hidden <- hidden + 1
          
        }
        
        
      }
      
      # To find how many seconds the snake does head-out after exiting the initial refuge
      if(sec_bf_body_out != time_limit){
        
        coordinates_body_out <- coordinates[sec_bf_body_out:nrow(coordinates), ]
        
        head_out_after_body_out <- 0
        
        for(j in 1:nrow(coordinates_body_out)){
          
          if(coordinates_body_out$behaviour[j] == "head-out"){
            
            head_out_after_body_out <- head_out_after_body_out + 1
            
          } else {
            
          }
          
        }
        
      } else {
        
        coordinates_body_out <- coordinates[1, ]
          
        head_out_after_body_out <- 0
        
      }

      
      
      
      tracking_ind$beh_moving <- moving
      
      tracking_ind$beh_exploring_wall <- exploring_wall
      
      tracking_ind$beh_immobile_exposed <- immobile
      
      tracking_ind$`beh_head-out` <- head_out
      
      tracking_ind$`beh_hiding-head-out` <- hidden_head_out
      
      tracking_ind$`beh_head-out_after_body_out` <- head_out_after_body_out
      
      tracking_ind$avg_speed <- mean(as.numeric(coordinates$dist_traveled))
      
      tracking_ind$total_dist_travelled <- sum(as.numeric(coordinates$dist_traveled))
      
      # Distance travelled when moving in the floor, not wall exploration
      moving <- coordinates[coordinates$behaviour == "moving", ]
      
      tracking_ind$avg_speed_exploring <- mean(as.numeric(moving$dist_traveled))
      
      tracking_ind$dist_travelled_exploring <- sum(as.numeric(moving$dist_traveled))
      
      # if the snake is out for more than 5 min
      if(nrow(coordinates_body_out) >= 300){
        
        # distance travelled in those 5 min
        tracking_ind$dist_trav_5min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:300]))
        
        # Then, if the snake has been out for 10 min or more
        if(nrow(coordinates_body_out) >= 600){
          
          # distance travelled in those 10 min
          tracking_ind$dist_trav_10min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:600]))
          
          # and if it has been outside for more than 15 min
          if(nrow(coordinates_body_out) >= 900){
            
            # distance travelled in those 15 min
            tracking_ind$dist_trav_15min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:900]))
            
            if(nrow(coordinates_body_out) >= 1200){
              
              # distance travelled in those 20 min
              tracking_ind$dist_trav_20min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:1200]))
              
              if(nrow(coordinates_body_out) >= 1500){
                
                # distance travelled in those 25 min
                tracking_ind$dist_trav_25min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:1500]))
                
                if(nrow(coordinates_body_out) >= 1800){
                  
                  # distance travelled in those 25 min
                  tracking_ind$dist_trav_30min_body_out <- sum(as.numeric(coordinates_body_out$dist_traveled[1:1800]))
                  
                }
                
                
              }
              
              
            }
            
          }
          
        }
        
      } else {
        
        # do nothing
        
      }
      
      
      # Merge ind data to global database
      
      summary_tracking[i, ] <- tracking_ind
      
      summary_tracking[i, ][summary_tracking[i, ] == "NaN"] <- NA
      
      summary_tracking[i, ][is.na(summary_tracking[i, ])] <- 0
      
      assign("summary_tracking", summary_tracking, envir=globalenv())
      
      
      
    } else {
      
      
      tracking_ind <- summary_tracking[1, ]
      
      tracking_ind[1, ] <- NA
      
      # ID
      tracking_ind$ID <- Snake_data_tracked$Snake_ID[i]
      
      # Area use
      tracking_ind$time_central_area <- 0
      
      tracking_ind$time_intermediate_area <- 0
      
      tracking_ind$time_refuge_area_1 <- 0
      
      tracking_ind$time_refuge_area_2 <- 0
      
      tracking_ind$time_initial_refuge <- 0
      
      tracking_ind$time_distant_wall <- 0
      
      tracking_ind$time_lateral_1 <- 0
      
      tracking_ind$time_lateral_2 <- 0
      
      tracking_ind$time_hidden_bf_head_out <- time_limit
      
      tracking_ind$total_time_hidden <- time_limit
      
      # Proximity to refuges
      tracking_ind$time_entrance <- 0
      
      tracking_ind$time_hole_1 <- 0
      
      tracking_ind$time_hole_2 <- 0
      
      tracking_ind$time_hole_3 <- 0
      
      tracking_ind$time_hole_4 <- 0
      
      tracking_ind$mean_dist_refuge <- 0
      
      tracking_ind$mean_dist_closest_refuge <- 0
      
      # Behaviour
      tracking_ind$beh_moving <- 0
      
      tracking_ind$beh_exploring_wall <- 0
      
      tracking_ind$beh_immobile_exposed <- 0
      
      tracking_ind$`beh_head-out` <- 0
      
      tracking_ind$`beh_hiding-head-out` <- 0
      
      tracking_ind$`beh_head-out_after_body_out` <- 0
      
      tracking_ind$avg_speed <- 0
      
      tracking_ind$avg_speed_exploring <- 0
      
      tracking_ind$total_dist_travelled <- 0
      
      tracking_ind$dist_travelled_exploring <- 0
      
      tracking_ind$dist_trav_5min_body_out <- 0
      
      tracking_ind$dist_trav_10min_body_out <- 0
      
      tracking_ind$dist_trav_15min_body_out <- 0
      
      tracking_ind$dist_trav_20min_body_out <- 0
      
      tracking_ind$dist_trav_25min_body_out <- 0
      
      tracking_ind$dist_trav_30min_body_out <- 0
      
      
      # Merge ind data to global database
      
      summary_tracking[i, ] <- tracking_ind
      
      assign("summary_tracking", summary_tracking, envir=globalenv())
      
    }
    
  }
  
}

generate_summary_tracking(3000)

# Save file in correspondant folder
write_xlsx(summary_tracking, "C:\\Users\\marc9\\Desktop\\Marc\\CREAF\\Snake tracking\\summary_tracking.xlsx")

#####


# Part 5: Merging data frames (already done)
#####

time_limit <- 3000

# Select columns of interest
Snake_data_tracked_final <- Snake_data_tracked[, c("Snake_ID", "Trap_ID", "Locality", "X", "Y", "Inv_situation", "Capture_day", "Test_day", "Acc_time", "Hour manip.", "Hour accl.", "Hour test", "Time_head_out_sec_raw", "Time_body_out_sec_raw", "Temperature", "R_H", "Weather", "Arena", "Official", "Repeatability", "Tracking", "Year_invasion", "Distance_Noahs", "SVL", "Real_Weight", "Condition", "Sex", "N_eggs", "Hard_Drive")]

# Transform variable into numeric
Snake_data_tracked_final$Time_head_out_sec_raw <- as.numeric(Snake_data_tracked_final$Time_head_out_sec_raw)

Snake_data_tracked_final$Time_body_out_sec_raw <- as.numeric(Snake_data_tracked_final$Time_body_out_sec_raw)

# If snake didn't get out, the number of seconds is 3000
Snake_data_tracked_final <- Snake_data_tracked_final %>% mutate(Time_head_out_sec_raw = ifelse(is.na(Time_head_out_sec_raw), time_limit, Time_head_out_sec_raw))

Snake_data_tracked_final <- Snake_data_tracked_final %>% mutate(Time_body_out_sec_raw = ifelse(is.na(Time_body_out_sec_raw), time_limit, Time_body_out_sec_raw))

# If head-out/Body-out time is larger than 3000 sec, cut it at 3000 
for(i in 1:nrow(Snake_data_tracked_final)){
  
  if(Snake_data_tracked_final$Time_head_out_sec_raw[i] > time_limit){
    
    Snake_data_tracked_final$Time_head_out_sec_raw[i] <- time_limit
    
  } else {
    
  }
  
  if(Snake_data_tracked_final$Time_body_out_sec_raw[i] > time_limit){
    
    Snake_data_tracked_final$Time_body_out_sec_raw[i] <- time_limit
    
  } else {
    
  }
  
}


# Change type of variable for some variables
Snake_data_tracked_final$Locality <- as.factor(Snake_data_tracked_final$Locality)
Snake_data_tracked_final$Inv_situation <- as.factor(Snake_data_tracked_final$Inv_situation)
Snake_data_tracked_final$Capture_day <- as.factor(Snake_data_tracked_final$Capture_day)
Snake_data_tracked_final$Test_day <- as.factor(Snake_data_tracked_final$Test_day)
Snake_data_tracked_final$Acc_time <- as.factor(Snake_data_tracked_final$Acc_time)
Snake_data_tracked_final$Arena <- as.factor(Snake_data_tracked_final$Arena)
Snake_data_tracked_final$Official <- as.factor(Snake_data_tracked_final$Official)
Snake_data_tracked_final$Repeatability <- as.factor(Snake_data_tracked_final$Repeatability)
Snake_data_tracked_final$Tracking <- as.factor(Snake_data_tracked_final$Tracking)


summary_tracking_complete <- cbind(Snake_data_tracked_final, summary_tracking)

summary_tracking_complete <- add_column(summary_tracking_complete, Visual_exp = summary_tracking_complete$Time_body_out_sec_raw - summary_tracking_complete$Time_head_out_sec_raw, .after = summary_tracking_complete$Time_body_out_sec_raw)

# Save file as CSV in order to plot it in QGIS and extract year of invasion, distance to Noah's garden, ...
write.csv2(summary_tracking_complete, "C:\\Users\\marc9\\Desktop\\Marc\\CREAF\\Snake tracking\\summary_tracking_complete.csv", row.names = FALSE)

write_xlsx(summary_tracking_complete, "C:\\Users\\marc9\\Desktop\\Marc\\CREAF\\Snake tracking\\summary_tracking_complete.xlsx")


#####


# Part 6: Plots
#####

time_limit <- 3000

# Which days did I test the snakes?
barplot(table(as.Date(summary_tracking_complete$Test_day)))

# At which hours did I test the snakes?
plot(summary_tracking_complete$`Hour test`)

# Some plots: These should be done without repeated indiviuals, only with new ones
summary_tracking_complete_no_rep <- summary_tracking_complete[summary_tracking_complete$Repeatability == "N",]

# Also, for some individuals, we are only going to use individuals that actually got out of the refuge (Body-out != 3000)
summary_tracking_complete_no_rep_outside <- summary_tracking_complete_no_rep[summary_tracking_complete_no_rep$Time_body_out_sec_raw != time_limit, ]

# Create one without gravid individuals?
summary_tracking_complete_no_rep_outside_noeggs <- summary_tracking_complete_no_rep_outside[summary_tracking_complete_no_rep_outside$N_eggs == 0, ]

# Part 6.1: Preliminary plots
#####

# Head-out ~ Tº
ggplot(summary_tracking_complete_no_rep, aes(x = Temperature, y = time_hidden_bf_head_out)) +
  geom_point() + 
  geom_smooth(method=lm, se = T) + 
  xlab("Temperature") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

# Head-out ~ Tº (group)
ggplot(summary_tracking_complete_no_rep, aes(x = Temperature, y = time_hidden_bf_head_out, color = Inv_situation)) + 
  geom_point() + 
  geom_smooth(method=lm, se = T) + 
  xlab("Temperature") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

# Head-out ~ Hour manipulation
ggplot(summary_tracking_complete_no_rep, aes(x = `Hour manip.`, y = time_hidden_bf_head_out)) +
  geom_point() + 
  geom_smooth(method=lm, se = T) + 
  xlab("Hour") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

# Head-out ~ depending on group
ggplot(aes(y = Time_body_out_sec_raw, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) +
  geom_boxplot(aes(fill = Inv_situation)) + 
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Time (s) hidden before head-out") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") +
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Time central area ~ depending on group
ggplot(aes(y = time_central_area, x = Inv_situation), data = summary_tracking_complete_no_rep) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Time (s) spent in central area (exposed)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") +
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Head-out behaviour ~ depending on group
ggplot(aes(y = `beh_head-out`, x = Inv_situation), data = summary_tracking_complete_no_rep) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Head-out behaviour (s)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Total distanced travelled ~ depending on group
ggplot(aes(y = total_dist_travelled, x = Inv_situation), data = summary_tracking_complete_no_rep) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Distanced explored (m)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Distanced travelled while exploring ~ depending on group
ggplot(aes(y = dist_travelled_exploring, x = Inv_situation), data = summary_tracking_complete_no_rep) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Distanced explored (m)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)


# Distanced travelled while exploring only if body-out ~ depending on group
ggplot(aes(y = dist_travelled_exploring, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Distanced explored (m)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Distanced travelled only if body-out ~ Year
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = dist_travelled_exploring)) +
  geom_point() + 
  geom_smooth(method=lm, se = FALSE) + 
  xlab("Temperature") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

cor.test(summary_tracking_complete_no_rep_outside$total_dist_travelled, summary_tracking_complete_no_rep_outside$Year_invasion)

# Time between head-out and body-out ~ depending on group (only for the ones that went out of the initial refuge)
ggplot(aes(y = Visual_exp, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Gap time between head-out and body-out (s)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Time head-out ~ Year invasion
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = Visual_exp)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("year of invasion") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

summary(lm(summary_tracking_complete_no_rep_outside$Visual_exp ~ summary_tracking_complete_no_rep_outside$Year_invasion))

# Time head-out ~ Year invasion by group
ggplot(summary_tracking_complete_no_rep, aes(x = Year_invasion, y = time_hidden_bf_head_out, color = Inv_situation)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("year of invasion") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()


# Distance travelled every 5 min
summary_tracking_complete_no_rep_outside <- summary_tracking_complete_no_rep_outside %>% mutate_at(c("dist_trav_5min_body_out", "dist_trav_10min_body_out", "dist_trav_15min_body_out", "dist_trav_20min_body_out", "dist_trav_25min_body_out", "dist_trav_30min_body_out"), ~na_if(., 0))

summary_tracking_complete_no_rep_outside$diff_10_5min <- summary_tracking_complete_no_rep_outside$dist_trav_10min_body_out - summary_tracking_complete_no_rep_outside$dist_trav_5min_body_out

summary_tracking_complete_no_rep_outside$diff_15_10min <- summary_tracking_complete_no_rep_outside$dist_trav_15min_body_out - summary_tracking_complete_no_rep_outside$dist_trav_10min_body_out

summary_tracking_complete_no_rep_outside$diff_20_15min <- summary_tracking_complete_no_rep_outside$dist_trav_20min_body_out - summary_tracking_complete_no_rep_outside$dist_trav_15min_body_out

summary_tracking_complete_no_rep_outside$diff_25_20min <- summary_tracking_complete_no_rep_outside$dist_trav_25min_body_out - summary_tracking_complete_no_rep_outside$dist_trav_20min_body_out

summary_tracking_complete_no_rep_outside$diff_30_25min <- summary_tracking_complete_no_rep_outside$dist_trav_30min_body_out - summary_tracking_complete_no_rep_outside$dist_trav_25min_body_out

distance_over_time <- tidyr::pivot_longer(summary_tracking_complete_no_rep_outside, 
                                       cols = starts_with("dist_trav_5min_body_out") | 
                                         starts_with("diff_10_5min") | 
                                         starts_with("diff_15_10min") | 
                                         starts_with("diff_20_15min") | 
                                         starts_with("diff_25_20min") | 
                                         starts_with("diff_30_25min"),
                                       names_to = "Time_Period",
                                       values_to = "Distance")

distance_over_time <- distance_over_time %>%
  mutate(Time_Period = factor(Time_Period, levels = c("dist_trav_5min_body_out", "diff_10_5min", "diff_15_10min", "diff_20_15min", "diff_25_20min", "diff_30_25min"),
                              labels = c("Δ 0-5 min", "Δ 5-10 min", "Δ 10-15 min", "Δ 15-20 min", "Δ 20-25 min", "Δ 25-30 min")))

ggplot(distance_over_time, aes(x = Time_Period, y = Distance, group = Snake_ID, color = Inv_situation)) +
  theme_bw() +
  geom_point(size = 1.5, position = position_nudge(x = 0)) +
  geom_line(size = 0.05, linetype = "dotted") +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, 
               position = position_nudge(x = 0.15)) + 
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", size = 0.2,
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.y = mean, geom = "line", size = 1, linetype = "solid", 
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.y = mean, geom = "point", size = 3, 
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, 
               position = position_nudge(x = -0.15)) + 
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", size = 0.2,
               position = position_nudge(x = -0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.y = mean, geom = "line", size = 1, linetype = "solid", 
               position = position_nudge(x = -0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.y = mean, geom = "point", size = 3, 
               position = position_nudge(x = -0.15)) +
  scale_x_discrete(drop = FALSE) +
  ylab("Distance (m)") +
  xlab("Time") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Distance traveled over time")


# Cumulative distance over time 
distance_over_time_sum <- tidyr::pivot_longer(summary_tracking_complete_no_rep_outside, 
                                          cols = starts_with("dist_trav_5min_body_out") | 
                                            starts_with("dist_trav_10min_body_out") | 
                                            starts_with("dist_trav_15min_body_out") | 
                                            starts_with("dist_trav_20min_body_out") | 
                                            starts_with("dist_trav_25min_body_out") | 
                                            starts_with("dist_trav_30min_body_out"),
                                          names_to = "Time_Period",
                                          values_to = "Distance")

distance_over_time_sum <- distance_over_time_sum %>%
  mutate(Time_Period = factor(Time_Period, levels = c("dist_trav_5min_body_out", "dist_trav_10min_body_out", "dist_trav_15min_body_out", "dist_trav_20min_body_out", "dist_trav_25min_body_out", "dist_trav_30min_body_out"),
                              labels = c("5 min", "10 min", "15 min", "20 min", "25 min", "30 min")))


ggplot(distance_over_time_sum, aes(x = Time_Period, y = Distance, group = Snake_ID, color = Inv_situation)) +
  theme_bw() +
  geom_point(size = 1.5, position = position_nudge(x = 0)) +
  geom_line(size = 0.05, linetype = "dotted") +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, 
               position = position_nudge(x = 0.15)) + 
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", size = 0.2,
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.y = mean, geom = "line", size = 1, linetype = "solid", 
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "CORE"), aes(group = Inv_situation),
               fun.y = mean, geom = "point", size = 3, 
               position = position_nudge(x = 0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, 
               position = position_nudge(x = -0.15)) + 
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", size = 0.2,
               position = position_nudge(x = -0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.y = mean, geom = "line", size = 1, linetype = "solid", 
               position = position_nudge(x = -0.15)) +
  stat_summary(data = . %>% filter(Inv_situation == "FRONT"), aes(group = Inv_situation),
               fun.y = mean, geom = "point", size = 3, 
               position = position_nudge(x = -0.15)) +
  scale_x_discrete(drop = FALSE) +
  ylab("Distance (m)") +
  xlab("Time") +
  ggtitle("Distance traveled over time")


##### 

# 6.1.1: Residual extraction
#####

residual_extraction_no_rep_outside <- function(){
  
  ## Condition
  
  # Head-out
  
  lm_head_out_cond <- lm(log1p(Time_head_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep_outside)
  
  summary_tracking_complete_no_rep_outside$res_head_cond <- lm_head_out_cond$residuals
  
  # by sex
  
  lm_head_out_cond_sex <- lm(log1p(Time_head_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside)
  
  summary_tracking_complete_no_rep_outside$res_head_cond_sex <- lm_head_out_cond_sex$residuals
  
  
  # Body-out
  
  lm_body_out_cond <- lm(log1p(Time_body_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_body_cond <- lm_body_out_cond$residuals
  
  # by sex
  
  lm_body_out_cond_sex <- lm(log1p(Time_body_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside)
  
  summary_tracking_complete_no_rep_outside$res_body_cond_sex <- lm_body_out_cond_sex$residuals
  
  
  # Visual exploration
  
  lm_visual_exp_cond <- lm(log1p(Visual_exp) ~ (Condition), data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_visual_cond <- lm_visual_exp_cond$residuals
  
  # lm by sex
  lm_visual_cond_sex <- lm(log1p(Visual_exp) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_visual_cond_sex <- lm_visual_cond_sex$residuals
  
  
  ## SVL
  
  # lm
  lm_head_out_svl <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_head_svl <- lm_head_out_svl$residuals
  
  # lm by sex
  lm_head_svl_sex <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_head_svl_sex <- lm_head_svl_sex$residuals
  
  
  # Body-out
  
  # lm
  lm_body_out_svl <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_body_svl <- lm_body_out_svl$residuals
  
  # lm by sex
  lm_body_svl_sex <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_body_svl_sex <- lm_body_svl_sex$residuals

  
  
  # Visual exploration
  
  lm_visual_exp_svl <- lm(log1p(Visual_exp) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_visual_svl <- lm_visual_exp_svl$residual
  
  # lm by sex
  lm_visual_svl_sex <- lm(log1p(Visual_exp) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside)

  summary_tracking_complete_no_rep_outside$res_visual_svl_sex <- lm_visual_svl_sex$residuals
  
  assign("summary_tracking_complete_no_rep_outside", summary_tracking_complete_no_rep_outside, envir=globalenv())

  
}

residual_extraction_no_rep_outside()

residual_extraction_no_rep <- function(){
  
  ## Condition
  
  # Head-out
  
  lm_head_out_cond <- lm(log1p(Time_head_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_head_cond <- lm_head_out_cond$residuals
  
  # by sex
  
  lm_head_out_cond_sex <- lm(log1p(Time_head_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_head_cond_sex <- lm_head_out_cond_sex$residuals
  
  
  # Body-out
  
  lm_body_out_cond <- lm(log1p(Time_body_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_body_cond <- lm_body_out_cond$residuals
  
  # by sex
  
  lm_body_out_cond_sex <- lm(log1p(Time_body_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_body_cond_sex <- lm_body_out_cond_sex$residuals
  
  
  # Visual exploration
  
  lm_visual_exp_cond <- lm(log1p(Visual_exp) ~ (Condition), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_visual_cond <- lm_visual_exp_cond$residuals
  
  # lm by sex
  lm_visual_cond_sex <- lm(log1p(Visual_exp) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_visual_cond_sex <- lm_visual_cond_sex$residuals
  
  
  ## SVL
  
  # lm
  lm_head_out_svl <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_head_svl <- lm_head_out_svl$residuals
  
  # lm by sex
  lm_head_svl_sex <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_head_svl_sex <- lm_head_svl_sex$residuals
  
  
  # Body-out
  
  # lm
  lm_body_out_svl <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_body_svl <- lm_body_out_svl$residuals
  
  # lm by sex
  lm_body_svl_sex <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_body_svl_sex <- lm_body_svl_sex$residuals
  
  
  
  # Visual exploration
  
  lm_visual_exp_svl <- lm(log1p(Visual_exp) ~ log1p(SVL), data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_visual_svl <- lm_visual_exp_svl$residual
  
  # lm by sex
  lm_visual_svl_sex <- lm(log1p(Visual_exp) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep)
  
  summary_tracking_complete_no_rep$res_visual_svl_sex <- lm_visual_svl_sex$residuals
  
  assign("summary_tracking_complete_no_rep", summary_tracking_complete_no_rep, envir=globalenv())
  
  
}

residual_extraction_no_rep()

residual_extraction_no_rep_outside_noeggs <- function(){
  
  ## Condition
  
  # Head-out
  
  lm_head_out_cond <- lm(log1p(Time_head_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_head_cond <- lm_head_out_cond$residuals
  
  # by sex
  
  lm_head_out_cond_sex <- lm(log1p(Time_head_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_head_cond_sex <- lm_head_out_cond_sex$residuals
  
  
  # Body-out
  
  lm_body_out_cond <- lm(log1p(Time_body_out_sec_raw) ~ (Condition), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_body_cond <- lm_body_out_cond$residuals
  
  # by sex
  
  lm_body_out_cond_sex <- lm(log1p(Time_body_out_sec_raw) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_body_cond_sex <- lm_body_out_cond_sex$residuals
  
  
  # Visual exploration
  
  lm_visual_exp_cond <- lm(log1p(Visual_exp) ~ (Condition), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_visual_cond <- lm_visual_exp_cond$residuals
  
  # lm by sex
  lm_visual_cond_sex <- lm(log1p(Visual_exp) ~ (Condition) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_visual_cond_sex <- lm_visual_cond_sex$residuals
  
  
  ## SVL
  
  # lm
  lm_head_out_svl <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_head_svl <- lm_head_out_svl$residuals
  
  # lm by sex
  lm_head_svl_sex <- lm(log1p(Time_head_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_head_svl_sex <- lm_head_svl_sex$residuals
  
  
  # Body-out
  
  # lm
  lm_body_out_svl <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_body_svl <- lm_body_out_svl$residuals
  
  # lm by sex
  lm_body_svl_sex <- lm(log1p(Time_body_out_sec_raw) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_body_svl_sex <- lm_body_svl_sex$residuals
  
  
  
  # Visual exploration
  
  lm_visual_exp_svl <- lm(log1p(Visual_exp) ~ log1p(SVL), data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_visual_svl <- lm_visual_exp_svl$residual
  
  # lm by sex
  lm_visual_svl_sex <- lm(log1p(Visual_exp) ~ log1p(SVL) + Sex, data = summary_tracking_complete_no_rep_outside_noeggs)
  
  summary_tracking_complete_no_rep_outside_noeggs$res_visual_svl_sex <- lm_visual_svl_sex$residuals
  
  assign("summary_tracking_complete_no_rep_outside_noeggs", summary_tracking_complete_no_rep_outside_noeggs, envir=globalenv())
  
  
}

residual_extraction_no_rep_outside_noeggs()

# binary sex

for(i in 1:nrow(summary_tracking_complete_no_rep_outside)){
  
  if(summary_tracking_complete_no_rep_outside$Sex[i] == "F"){
    
    summary_tracking_complete_no_rep_outside$Sex_12[i] <- 1
    
  } else if(summary_tracking_complete_no_rep_outside$Sex[i] == "M"){
    
    summary_tracking_complete_no_rep_outside$Sex_12[i] <- 2
    
  } else {
    
    summary_tracking_complete_no_rep_outside$Sex_12[i] <- 3
    
  }
  
}


# plots, together and by sex
ggplot(summary_tracking_complete_no_rep_outside, aes(x = (Condition), y = log1p(Time_head_out_sec_raw))) +
  geom_point() +
  stat_smooth(method = "lm") + 
  theme_minimal()

ggplot(summary_tracking_complete_no_rep_outside, aes(x = (Condition), y = log1p(Time_head_out_sec_raw), color = Sex)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  theme_minimal()



ggplot(summary_tracking_complete_no_rep_outside, aes(x = Distance_Noahs, y = res_visual_svl_sex, col = Sex)) +
  geom_point() + 
  geom_smooth(method=lm, se = T) + 
  xlab("Distance to Noah's Garden") + 
  ylab("Body-out residuals") + 
  theme_minimal()



summary(lm(res_visual_svl_sex ~ Year_invasion * Sex, data = summary_tracking_complete_no_rep_outside))


summary(lm(res_visual_svl_sex ~ Year_invasion, data = summary_tracking_complete_no_rep_outside, subset = (Sex_12 == 1)))


#####

# Part 6.2: Survival analysis
#####

# Survival analysis for head-out and body-out analysis

# We add columns for saying if data is censored (0) or the event has happened within the time frame (1)
summary_tracking_complete_no_rep <- add_column(summary_tracking_complete_no_rep, Head_out_01 = NA, .after = "Time_head_out_sec_raw")


summary_tracking_complete_no_rep <- add_column(summary_tracking_complete_no_rep, Body_out_01 = NA, .after = "Time_body_out_sec_raw")


# 1 for event (head-out/body-out) and 0 for censored data
for(i in 1:nrow(summary_tracking_complete_no_rep)){
  
  if(summary_tracking_complete_no_rep$Time_head_out_sec_raw[i] == 3000){
    
    summary_tracking_complete_no_rep$Head_out_01[i] <- 0
    
  } else {
    
    summary_tracking_complete_no_rep$Head_out_01[i] <- 1
    
  }
  
  if(summary_tracking_complete_no_rep$Time_body_out_sec_raw[i] == 3000){
    
    summary_tracking_complete_no_rep$Body_out_01[i] <- 0
    
  } else {
    
    summary_tracking_complete_no_rep$Body_out_01[i] <- 1
    
  }
  
}

# + indicates censored data
Surv(summary_tracking_complete_no_rep$Time_head_out_sec_raw, summary_tracking_complete_no_rep$Head_out_01)


survival_curve_1 <- survfit(Surv(Time_head_out_sec_raw, Head_out_01) ~ 1,
                            data = summary_tracking_complete_no_rep)

str(survival_curve_1)

plot(survival_curve_1)

survfit2(Surv(Time_head_out_sec_raw, Head_out_01) ~ 1,
         data = summary_tracking_complete_no_rep) %>% 
  ggsurvfit() +
  labs(
    x = "Time (s)",
    y = "Overall head-out probability"
  ) +
  add_confidence_interval() +
  add_risktable()

# Probability head-out within 30 min = 1 - 0.318 = 0.682
summary(survfit(Surv(Time_head_out_sec_raw, Head_out_01) ~ 1,
                data = summary_tracking_complete_no_rep), times = 1800)

# Table that gives us the probability of an event depending on if we use censored data or not
survfit(Surv(Time_head_out_sec_raw, Head_out_01) ~ 1,
        data = summary_tracking_complete_no_rep) %>% 
  tbl_survfit(
    times = 1800,
    label_header = "**30 min head-out (95% CI)**"
  )

# Head-out probability between front and core
survfit2(Surv(Time_head_out_sec_raw, Head_out_01) ~ Inv_situation,
         data = summary_tracking_complete_no_rep) %>% 
  ggsurvfit() +
  labs(
    x = "Time (s)",
    y = "Overall head-out probability"
  ) +
  add_confidence_interval() +
  add_risktable()

fit_head_out <- survfit(Surv(Time_head_out_sec_raw, Head_out_01) ~ Inv_situation, data = summary_tracking_complete_no_rep)


ggsurvplot(fit_head_out,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = c(1, 1),
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           surv.median.line = "hv",
           legend.labs = c("CORE", "FRONT"))


# Body-out probability between front and core
survfit2(Surv(Time_body_out_sec_raw, Body_out_01) ~ Inv_situation,
         data = summary_tracking_complete_no_rep) %>% 
  ggsurvfit() +
  labs(
    x = "Time (s)",
    y = "Overall Body-out probability"
  ) +
  add_confidence_interval() +
  add_risktable()


fit_body_out <- survfit(Surv(Time_body_out_sec_raw, Body_out_01) ~ Inv_situation, data = summary_tracking_complete_no_rep)


ggsurvplot(fit_body_out,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = c(1, 1),
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           surv.median.line = "hv",
           legend.labs = c("CORE", "FRONT"))



# Smooth survival plot (general)
sm.options(
  list(
    xlab = "Invasion year",
    ylab = "Head-out")
)

sm.survival(
  x = summary_tracking_complete_no_rep$Year_invasion,
  y = summary_tracking_complete_no_rep$Time_head_out_sec_raw,
  status = summary_tracking_complete_no_rep$Head_out_01,
  h = sd(summary_tracking_complete_no_rep$Year_invasion) / nrow(summary_tracking_complete_no_rep)^(-1/4)
)

sm.options(
  list(
    xlab = "Invasion year",
    ylab = "Body-out")
)

sm.survival(
  x = summary_tracking_complete_no_rep$Year_invasion,
  y = summary_tracking_complete_no_rep$Time_body_out_sec_raw,
  status = summary_tracking_complete_no_rep$Body_out_01,
  h = sd(summary_tracking_complete_no_rep$Year_invasion) / nrow(summary_tracking_complete_no_rep)^(-1/4)
)



#####


# Part 6.3: Repeatability analysis

# Part 6.3.1: Head-out
#####

# Repeatability

# Select rows with duplicated ID (repeatability)

repeat_ID <- summary_tracking_complete[duplicated(summary_tracking_complete$ID), "ID"]

# Select ind with more than one tracking
repeat_ind <- summary_tracking_complete[summary_tracking_complete$Snake_ID %in% repeat_ID, ]


# assign 1 to first trials
repeat_ind$repeat_trial <- 1

for (i in 1:nrow(repeat_ind)){
  
  # And if duplicated, assign 2
  if(duplicated(repeat_ind$ID)[i] == T){
    
    repeat_ind$repeat_trial[i] <- 2
    
  } else {
    
  }
  
}

# Order snakes by ID (first the first trial and second the repeatiblity test)
repeat_ind <- repeat_ind[order(repeat_ind$ID), ]

# empty plot 
plot(NULL, xlim = c(1, 2), ylim = c(0, 3000), xlab = "First trial (left)   |   Repeatibility trial (right)", ylab = "Time to head-out", xaxt = "n")

# connect pairs of dots (same ind, first and second trial)
for(i in 1:(nrow(repeat_ind)/2)){
  
  if(repeat_ind$Inv_situation[2*i] == "FRONT"){
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$time_hidden_bf_head_out[(2*i) - 1], repeat_ind$time_hidden_bf_head_out[(2*i)]), type = "b", col = "red")
    
  } else {
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$time_hidden_bf_head_out[(2*i) - 1], repeat_ind$time_hidden_bf_head_out[(2*i)]), type = "b", col = "blue")
    
  }
  
}

# Complete plot ~ invasion status
ggplot(data = repeat_ind, aes(x = Repeatability, y = Time_head_out_sec_raw, color = Inv_situation)) +
  theme_bw() + 
  stat_summary(
    geom = "errorbar", 
    fun.data = mean_cl_boot, 
    width = 0.1,
    size = 0.8,  
    col = "grey57") + 
  geom_point(size = 3) +
  stat_summary(fun.y = mean, position = "dodge", size = 1) +
  stat_summary(aes(group = Inv_situation), geom = "line", fun.y = mean, size = 1.5) +  
  geom_line(aes(group = Snake_ID), size = 0.5) + 
  ylab("Time (s)") + 
  ggtitle("Head-out time (s) since start of the experiment")

# General plot (no groups)
ggplot(data = repeat_ind, aes(x = Repeatability, y = Time_head_out_sec_raw)) +
  theme_bw() + 
  stat_summary(
    geom = "errorbar", 
    fun.data = mean_cl_boot, 
    width = 0.1,
    size = 0.8,  
    col = "grey57") + 
  geom_point(size = 3) +
  stat_summary(fun.y = mean, position = "dodge", size = 1) +
  stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.5) +  
  geom_line(aes(group = Snake_ID), size = 0.5) + 
  ylab("Time (s)") + 
  ggtitle("Head-out time (s) since start of the experiment")


# interclass correlation indexes (ICC, package IRR)

# We extract the variable that we want to test repeatability
icc_data <- as.data.frame(repeat_ind[, c("Snake_ID", "time_hidden_bf_head_out", "repeat_trial")])

# Put it in wide format
icc_data_wide <- reshape(icc_data, idvar = "Snake_ID", timevar = "repeat_trial", direction = "wide")

# Delete ID column
icc_data_wide <- icc_data_wide[, -1]

# Change column names
colnames(icc_data_wide) <- c("First_trial", "Second_trial")

# Consistency (ICC)
icc(icc_data_wide, model = "oneway", type = "consistency", unit = "single")

# An intraclass correlation coefficient, according to Koo & Li:
# Less than 0.50: Poor reliability
# Between 0.5 and 0.75: Moderate reliability
# Between 0.75 and 0.9: Good reliability
# Greater than 0.9: Excellent reliability

# Alternative plot
# linear regression
summary(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))

# plot
plot(icc_data_wide$Second_trial ~ icc_data_wide$First_trial, xlab = "First trial's head-out time (s)", ylab = "Second trial's head-out time (s)")

abline(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))


#####


# Part 6.3.2: Body-out
#####

# empty plot 
plot(NULL, xlim = c(1, 2), ylim = c(0, 3000), xlab = "First trial (left)   |   Repeatibility trial (right)", ylab = "Time to body-out", xaxt = "n")

# connect pairs of dots (same ind, first and second trial)
for(i in 1:(nrow(repeat_ind)/2)){
  
  if(repeat_ind$Inv_situation[2*i] == "FRONT"){
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$Time_body_out_sec_raw[(2*i) - 1], repeat_ind$Time_body_out_sec_raw[(2*i)]), type = "b", col = "red")
    
  } else {
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$Time_body_out_sec_raw[(2*i) - 1], repeat_ind$Time_body_out_sec_raw[(2*i)]), type = "b", col = "blue")
    
  }
  
}


# interclass correlation indexes (ICC, package IRR)

# We extract the variable that we want to test repeatability
icc_data <- as.data.frame(repeat_ind[, c("Snake_ID", "Time_body_out_sec_raw", "repeat_trial")])

# Put it in wide format
icc_data_wide <- reshape(icc_data, idvar = "Snake_ID", timevar = "repeat_trial", direction = "wide")

# Delete ID column
icc_data_wide <- icc_data_wide[, -1]

# Change column names
colnames(icc_data_wide) <- c("First_trial", "Second_trial")

# Consistency (ICC)
icc(icc_data_wide, model = "oneway", type = "consistency", unit = "single")

# An intraclass correlation coefficient, according to Koo & Li:
# Less than 0.50: Poor reliability
# Between 0.5 and 0.75: Moderate reliability
# Between 0.75 and 0.9: Good reliability
# Greater than 0.9: Excellent reliability

# Alternative plot
# linear regression
summary(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))

# plot
plot(icc_data_wide$Second_trial ~ icc_data_wide$First_trial, xlab = "First trial's body-out time (s)", ylab = "Second trial's body-out time (s)")

abline(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))


#####


# Part 6.3.3: Visual exploration
#####

# empty plot 
plot(NULL, xlim = c(1, 2), ylim = c(0, 3000), xlab = "First trial (left)   |   Repeatibility trial (right)", ylab = "Gap time", xaxt = "n")

repeat_ind$Visual_exp <- repeat_ind$Time_body_out_sec_raw - repeat_ind$Time_head_out_sec_raw

# connect pairs of dots (same ind, first and second trial)
for(i in 1:(nrow(repeat_ind)/2)){
  
  if(repeat_ind$Inv_situation[2*i] == "FRONT"){
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$Visual_exp[(2*i) - 1], repeat_ind$Visual_exp[(2*i)]), type = "b", col = "red")
    
  } else {
    
    lines(c(repeat_ind$repeat_trial[(2*i) - 1], repeat_ind$repeat_trial[(2*i)]), c(repeat_ind$Visual_exp[(2*i) - 1], repeat_ind$Visual_exp[(2*i)]), type = "b", col = "blue")
    
  }
  
}


# interclass correlation indexes (ICC, package IRR)

# We extract the variable that we want to test repeatability
icc_data <- as.data.frame(repeat_ind[, c("Snake_ID", "Visual_exp", "repeat_trial")])

# Put it in wide format
icc_data_wide <- reshape(icc_data, idvar = "Snake_ID", timevar = "repeat_trial", direction = "wide")

# Delete ID column
icc_data_wide <- icc_data_wide[, -1]

# Change column names
colnames(icc_data_wide) <- c("First_trial", "Second_trial")

# Consistency (ICC)
icc(icc_data_wide, model = "oneway", type = "consistency", unit = "single")

# An intraclass correlation coefficient, according to Koo & Li:
# Less than 0.50: Poor reliability
# Between 0.5 and 0.75: Moderate reliability
# Between 0.75 and 0.9: Good reliability
# Greater than 0.9: Excellent reliability

# Alternative plot
# linear regression
summary(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))

# plot
plot(icc_data_wide$Second_trial ~ icc_data_wide$First_trial, xlab = "First trial's gap time (s)", ylab = "Second trial's gap time (s)")

abline(lm(icc_data_wide$Second_trial ~ icc_data_wide$First_trial))

#####









#####


# Part 6.3.4: Joint plot
#####

head_out_rep <- ggplot(data = repeat_ind, aes(x = Repeatability, y = Time_head_out_sec_raw)) +
  theme_bw() + 
  stat_summary(
    geom = "errorbar", 
    fun.data = mean_cl_boot, 
    width = 0.1,
    size = 0.8,  
    col = "black",
    position = position_nudge(x = c(-0.2, 0.2))) +
  geom_point(size = 3) +
  stat_summary(fun.y = mean, position = position_nudge(x = c(-0.2, 0.2)), size = 1) +
  stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.5, linetype = "solid", position = position_nudge(x = c(-0.2, 0.2))) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", color = "black", size = 0.2,
               position = position_nudge(x = c(-0.2, 0.2))) +
  geom_line(aes(group = Snake_ID), size = 0.5, linetype = "dashed") + 
  ylab("Time (s)") + 
  ggtitle("Head-out") + 
  labs(caption = "ICC = 0.855 | p < 0.001") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, NA))


body_out_rep <- ggplot(data = repeat_ind, aes(x = Repeatability, y = Time_body_out_sec_raw)) +
  theme_bw() + 
  stat_summary(
    geom = "errorbar", 
    fun.data = mean_cl_boot, 
    width = 0.1,
    size = 0.8,  
    col = "black",
    position = position_nudge(x = c(-0.2, 0.2))) +
  geom_point(size = 3) +
  stat_summary(fun.y = mean, position = position_nudge(x = c(-0.2, 0.2)), size = 1) +
  stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.5, linetype = "solid", position = position_nudge(x = c(-0.2, 0.2))) + 
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", color = "black", size = 0.2,
               position = position_nudge(x = c(-0.2, 0.2))) +
  geom_line(aes(group = Snake_ID), size = 0.5, linetype = "dashed") + 
  ylab("Time (s)") + 
  ggtitle("Body-out") + 
  labs(caption = "ICC = 0.777 | p < 0.001") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, NA))

repeat_ind$Visual_exp <- repeat_ind$Time_body_out_sec_raw - repeat_ind$Time_head_out_sec_raw

Visual_exp_rep <- ggplot(data = repeat_ind, aes(x = Repeatability, y = Visual_exp)) +
  theme_bw() + 
  stat_summary(
    geom = "errorbar", 
    fun.data = mean_cl_boot, 
    width = 0.1,
    size = 0.8,  
    col = "black",
    position = position_nudge(x = c(-0.2, 0.2))) + 
  geom_point(size = 3) +
  stat_summary(fun.y = mean, position = position_nudge(x = c(-0.2, 0.2)), size = 1) +
  stat_summary(aes(group = 1), fun.y = mean, geom = "line", size = 1.5, linetype = "solid", position = position_nudge(x = c(-0.2, 0.2))) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "pointrange", color = "black", size = 0.2,
               position = position_nudge(x = c(-0.2, 0.2))) +
  geom_line(aes(group = Snake_ID), size = 0.5, linetype = "dashed") +
  ylab("Time (s)") + 
  ggtitle("Visual exploration") + 
  labs(caption = "ICC = 0.382 | p = 0.0497") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, NA))


ggarrange(head_out_rep, body_out_rep, Visual_exp_rep, 
          ncol = 3, nrow = 1)


#####


# Part 6.4: Behaviour and area use
#####

# Behaviour
# Select columns of interest for stacked behaviour plot
stacked_behaviour <- summary_tracking_complete_no_rep[, c("Snake_ID", "time_hidden_bf_head_out", "total_time_hidden", "beh_moving", "beh_exploring_wall", "beh_immobile_exposed", "beh_head-out", "beh_hiding-head-out")]

# New variable: hidden time after first head-out
stacked_behaviour$`beh_hidden_after_first_head-out` <- stacked_behaviour$total_time_hidden - stacked_behaviour$time_hidden_bf_head_out

# We delete this intermediate variable used to calculate previous variable
stacked_behaviour <- subset(stacked_behaviour, select = -total_time_hidden)

# Change column names
colnames(stacked_behaviour) <- c("Snake_ID", "Hidden before head-out", "Exploring", "Exploring wall", "Immobile exposed", "Head-out", "Hidden in corner", "Hidden after first head-out")

# Alternative data frame with invasion situation
stacked_behaviour_treatment <- stacked_behaviour

stacked_behaviour_treatment$treatment <- summary_tracking_complete_no_rep$Inv_situation

stacked_behaviour_treatment$year <- summary_tracking_complete_no_rep$Year_invasion

# Snake_ID ordered by year of invasion
snakes_order_by_year <- summary_tracking_complete_no_rep[, c("Snake_ID", "Year_invasion")]

snakes_order_by_year <- snakes_order_by_year[order(snakes_order_by_year$Year_invasion, decreasing = FALSE), ]

# Long format
stacked_behaviour_long <- pivot_longer(stacked_behaviour, cols = c("Hidden before head-out", "Exploring", "Exploring wall", "Immobile exposed", "Head-out", "Hidden in corner", "Hidden after first head-out"), names_to = "Behaviour", values_to = "Percentage")

# Percentage of each behaviour
stacked_behaviour_long <- stacked_behaviour_long %>%
  group_by(Snake_ID) %>%
  mutate(Percentage_total = sum(Percentage)) %>%
  ungroup() %>%
  mutate(Percentage = Percentage / Percentage_total * 100)

# Order Snake ID by year of invasion
stacked_behaviour_long$Snake_ID <- factor(stacked_behaviour_long$Snake_ID, levels = snakes_order_by_year$Snake_ID)

# Order behaviour variables
stacked_behaviour_long$Behaviour <- factor(stacked_behaviour_long$Behaviour, levels = c("Immobile exposed", "Exploring wall", "Exploring", "Head-out", "Hidden in corner", "Hidden after first head-out", "Hidden before head-out"))

# Colors for behaviours
beh_colors <- c(viridis(7)[7], viridis(7)[6], viridis(7)[5], viridis(7)[4], viridis(7)[3], viridis(7)[2], viridis(7)[1])

# Plot
ggplot(stacked_behaviour_long, aes(x = Snake_ID, y = Percentage, fill = Behaviour)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Snake ID", y = "% behaviours during 50 min test", fill = "Behaviours") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
        legend.text = element_text(size = 10), 
        axis.title.x = element_text(size = 15),  
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 10), 
        legend.title = element_text(size = 15)) +
  scale_fill_manual(values = beh_colors)


# Area
#Select columns
stacked_area <- summary_tracking_complete_no_rep[, c("Snake_ID", "time_central_area", "time_intermediate_area", "time_refuge_area_1", "time_refuge_area_2", "time_initial_refuge", "time_distant_wall", "time_lateral_1", "time_lateral_2", "time_hidden_bf_head_out", "total_time_hidden")]

# New variables
stacked_area$refuge_area <- stacked_area$time_refuge_area_1 + stacked_area$time_refuge_area_2

stacked_area$over_lateral <- stacked_area$time_lateral_1 + stacked_area$time_lateral_2

stacked_area$`beh_hidden_after_first_head-out` <- stacked_area$total_time_hidden - stacked_area$time_hidden_bf_head_out

# We delete theese intermediate variables used to calculate previous variables
stacked_area <- subset(stacked_area, select = -c(total_time_hidden, time_refuge_area_1, time_refuge_area_2, time_lateral_1, time_lateral_2))

# Change column names
colnames(stacked_area) <- c("Snake_ID", "Central area", "Intermediate area", "Initial refuge", "Distant wall", "Hidden before head-out", "Refuge area", "Over lateral refuge", "Head-out after head-out")

# Alternative data frame with invasion situation
stacked_area_treatment <- stacked_area

stacked_area_treatment$treatment <- summary_tracking_complete_no_rep$Inv_situation

stacked_area_treatment$year <- summary_tracking_complete_no_rep$Year_invasion

# Long format
stacked_area_long <- pivot_longer(stacked_area, cols = c("Central area", "Intermediate area", "Initial refuge", "Distant wall", "Hidden before head-out", "Refuge area", "Over lateral refuge", "Head-out after head-out"), names_to = "Area", values_to = "Percentage")

# Percentage of each behaviour
stacked_area_long <- stacked_area_long %>%
  group_by(Snake_ID) %>%
  mutate(Percentage_total = sum(Percentage)) %>%
  ungroup() %>%
  mutate(Percentage = Percentage / Percentage_total * 100)

# Order Snake ID by year of invasion
stacked_area_long$Snake_ID <- factor(stacked_area_long$Snake_ID, levels = snakes_order_by_year$Snake_ID)

# Order behaviour variables
stacked_area_long$Area <- factor(stacked_area_long$Area, levels = c("Central area", "Intermediate area", "Initial refuge", "Distant wall", "Refuge area", "Over lateral refuge", "Head-out after head-out", "Hidden before head-out"))


# Colors for behaviours
area_colors <- c(viridis(8)[8], viridis(8)[7], viridis(8)[6], viridis(8)[5], viridis(8)[4], viridis(8)[3], viridis(8)[2], viridis(8)[1])

# Plot
ggplot(stacked_area_long, aes(x = Snake_ID, y = Percentage, fill = Area)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Snake ID", y = "% area during 50 min test", fill = "Area") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
        legend.text = element_text(size = 10), 
        axis.title.x = element_text(size = 15),  
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 10), 
        legend.title = element_text(size = 15)) +
  scale_fill_manual(values = area_colors)


# Spatial analysis + PERMANOVA / SIMPER (ANOVA?)

# Behaviour
# Delete Snake ID column
stacked_behaviour_treatment <- subset(stacked_behaviour_treatment, select = -c(Snake_ID))

# Subset pca columns (numeric)
stacked_behaviour_treatment_var <- stacked_behaviour_treatment[1:7]

# PCA
pca_beh <- prcomp(stacked_behaviour_treatment_var, scale. = T)

# Weights PCA axis
pca_beh$rotation

# Coordinates
pca_beh$x

# Summary
summary(pca_beh)

# Plot 1
autoplot(pca_beh, data = stacked_behaviour_treatment, colour = 'treatment',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# Plot 2
fviz_pca_biplot(pca_beh, geom.ind = "point", 
             habillage = stacked_behaviour_treatment$treatment, 
             axes = c(1, 2), 
             pointsize = 3,
             addEllipses = TRUE, ellipse.level = 0.9) 


# Area

# Delete Snake ID column
stacked_area_treatment <- subset(stacked_area_treatment, select = -c(Snake_ID))

# Subset pca columns (numeric)
stacked_area_treatment_var <- stacked_area_treatment[, c("Central area", "Intermediate area", "Initial refuge", "Distant wall", "Refuge area", "Over lateral refuge")]

# PCA
pca_area <- prcomp(stacked_area_treatment_var, scale. = T)

# Weights PCA axis
pca_area$rotation

# Coordinates
pca_area$x

# Summary
summary(pca_area)

# Plot 1
autoplot(pca_area, data = stacked_area_treatment, colour = 'treatment',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# Plot 2
fviz_pca_biplot(pca_area, geom.ind = "point", 
                habillage = stacked_area_treatment$treatment, 
                axes = c(1, 2), 
                pointsize = 3,
                addEllipses = TRUE, ellipse.level = 0.9) 


# PERMANOVA
stacked_area_treatment_env <- stacked_area_treatment[, c("treatment", "year")]

stacked_behaviour_treatment_env <- stacked_behaviour_treatment[, c("treatment", "year")]

# PERMANOVA behaviour
beh.permanova <- adonis2(stacked_behaviour_treatment_var ~ treatment + year, data = stacked_behaviour_treatment_env, permutations = 999, method="euclidean")

beh.permanova

# PERMANOVA area
beh.area <- adonis2(stacked_area_treatment_var ~ treatment + year, data = stacked_area_treatment_env, permutations = 999, method="euclidean")

beh.area


# Alternative method behaviour
beh.dist <- vegdist(stacked_behaviour_treatment_var, method="euclidean")

beh.dist <- vegdist(decostand(stacked_behaviour_treatment_var, "norm"), method="euclidean")


beh.dispersion <- betadisper(beh.dist, group = stacked_behaviour_treatment_env$treatment)

permutest(beh.dispersion)

plot(beh.dispersion, hull = FALSE, ellipse = TRUE) 



# Alternative method area
area.dist <- vegdist(stacked_area_treatment_var, method="euclidean")

area.dist <- vegdist(decostand(stacked_area_treatment_var, "norm"), method="euclidean")

area.dispersion <- betadisper(area.dist, group = stacked_area_treatment_env$treatment)

permutest(area.dispersion)

plot(area.dispersion, hull = FALSE, ellipse = TRUE) 




 #####



# Part 7: Plots by invasion category
#####

# Head-out
ggsurvplot(fit_head_out,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = c(1, 1),
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           surv.median.line = "hv",
           legend.labs = c("CORE", "FRONT"))


# Alternative method area
area.dist <- vegdist(stacked_area_treatment_var, method="euclidean")

beh.dist <- vegdist(decostand(stacked_area_treatment_var, "norm"), method="euclidean")


beh.dispersion <- betadisper(beh.dist, group = stacked_area_treatment_env$treatment)

permutest(beh.dispersion)

plot(beh.dispersion, hull = FALSE, ellipse = TRUE) ggplot(aes(y = time_hidden_bf_head_out, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) +
  geom_boxplot(aes(fill = Inv_situation)) + 
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Time (s) hidden before head-out") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") +
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)


# Body-out
ggsurvplot(fit_body_out,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = c(1, 1),
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           surv.median.line = "hv",
           legend.labs = c("CORE", "FRONT"))


ggplot(aes(y = Time_body_out_sec_raw , x = Inv_situation), data = summary_tracking_complete_no_rep_outside) +
  geom_boxplot(aes(fill = Inv_situation)) + 
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Time (s) hidden before body-out") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") +
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)



# Visual exploration
ggplot(aes(y = Time_body_out_sec_raw - Time_head_out_sec_raw , x = Inv_situation), data = summary_tracking_complete_no_rep_outside) +
  geom_boxplot(aes(fill = Inv_situation)) + 
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Time (s) doing visual exploration") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") +
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Distance explored
ggplot(aes(y = dist_travelled_exploring, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Distanced explored (m)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)


# Head-out after body out
ggplot(aes(y = `beh_head-out_after_body_out`, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Head-out (s) after first emergence") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)


# Time hidden after first head-out
ggplot(aes(y = total_time_hidden - time_hidden_bf_head_out, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Head-out (s) after first emergence") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Time in central area
ggplot(aes(y = time_central_area + time_intermediate_area, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Central area (s)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

# Time exploring
ggplot(aes(y = beh_moving, x = Inv_situation), data = summary_tracking_complete_no_rep_outside) + 
  geom_boxplot(aes(fill = Inv_situation)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) + 
  xlab("Invasion status") + 
  ylab("Exploring (s)") + 
  stat_summary(fun=mean, geom="point", shape=21, size=4, color="black", fill="white") + 
  theme_minimal() + 
  stat_n_text(y.pos = NULL, color = "black", text.box = TRUE)

#####


# Part 8: Linear plots 
#####

# Head-out
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = time_hidden_bf_head_out)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Time (s) hidden before head-out") + 
  theme_minimal()

# Body-out
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = Time_body_out_sec_raw)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Time (s) hidden before body-out") + 
  theme_minimal()

# Visual exploration
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = Time_body_out_sec_raw - Time_head_out_sec_raw)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Visual explortion") + 
  theme_minimal()

# Distance explored
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = dist_travelled_exploring)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Distance explored") + 
  theme_minimal()

# Exploring
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = beh_moving)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Behaviour moving") + 
  theme_minimal()

# Immobile exposed
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = beh_immobile_exposed)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Behaviour exposed") + 
  theme_minimal()

# Head-out after head-out
ggplot(summary_tracking_complete_no_rep_outside, aes(x = Year_invasion, y = `beh_head-out_after_body_out`)) +
  geom_point() + 
  geom_smooth(method=lm) + 
  xlab("Year of invasion") + 
  ylab("Head-out after body out") + 
  theme_minimal()




#####



##### Part 9: Maximum likelihood estimates regression (to deal with censored data) #####

# Example
# https://www.cfholbert.com/blog/censored-regression/

install.packages("magrittr")
install.packages("censReg")
install.packages("summarytools")

library(magrittr)
library(censReg)
library(summarytools)


# Head-out

# Plot censored data (only right cenosred, maximum time it takes the snake to exit the refuge)
plot(Time_head_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep,  pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')

# Censored regression
fit.cen <- censReg(Time_head_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep, 
                   left = -Inf, right = 3000)
summary(fit.cen)
# There is no significant regression (p = 0.505)

# Get the standard error. It is pretty big (more than a thousand)
coef(fit.cen, logSigma = FALSE)


# Fit a line using OLS regression (linear regression)
fit.OLS <- lm(Time_head_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep)
summary(fit.OLS)

# fit regression only of non-censored data
fit.detect <- lm(summary_tracking_complete_no_rep$Time_head_out_sec_raw[0 <= summary_tracking_complete_no_rep$Time_head_out_sec_raw & summary_tracking_complete_no_rep$Time_head_out_sec_raw < 3000] ~ summary_tracking_complete_no_rep$Year_invasion[0 <= summary_tracking_complete_no_rep$Time_head_out_sec_raw & summary_tracking_complete_no_rep$Time_head_out_sec_raw < 3000])
summary(fit.detect)

# Table with slope estimates of each linear fit
rbind(
  Censored = c(coef(fit.cen)[1], coef(fit.cen)[2]),
  OLS_Full = coef(fit.OLS),
  OLS_Detects = coef(fit.detect)
) %>% round(4)

# Plot the results
lineplot <- function() {
  
  # Add lines
  abline(coef(fit.cen)[1:2], col = 'Red', lwd = 2)
  abline(coef(fit.OLS), col = 'Blue', lty = 2, lwd = 2)
  abline(coef(fit.detect), col = rgb(0.2, 0.6, 0.2), lty = 3, lwd = 2)
  
  # Add legend
  legend(
    "topleft",
    legend = c(
      "Censored",
      "OLS",
      "OLS Detects"
    ),
    lwd = 3 * par("cex"),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 3),
    text.font = 1, bty = "n",
    pt.cex = 1, cex = 0.8, y.intersp = 1.3
  )
}

par(mfrow = c(1, 2))

plot(summary_tracking_complete_no_rep$Year_invasion, summary_tracking_complete_no_rep$Time_head_out_sec_raw, pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')
lineplot()

hist(summary_tracking_complete_no_rep$Time_head_out_sec_raw, breaks = 100, main = 'Censored Data')



# Body out

plot(Time_body_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep,  pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')

# Censored regression
fit.cen <- censReg(Time_body_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep, 
                   left = -Inf, right = 3000)
summary(fit.cen)


# Get the standard error. It is pretty big (more than a thousand)
coef(fit.cen, logSigma = FALSE)


# Fit a line using OLS regression (linear regression)
fit.OLS <- lm(Time_body_out_sec_raw ~ Year_invasion, data = summary_tracking_complete_no_rep)
summary(fit.OLS)

# fit regression only of non-censored data
fit.detect <- lm(summary_tracking_complete_no_rep$Time_body_out_sec_raw[0 <= summary_tracking_complete_no_rep$Time_body_out_sec_raw & summary_tracking_complete_no_rep$Time_body_out_sec_raw < 3000] ~ summary_tracking_complete_no_rep$Year_invasion[0 <= summary_tracking_complete_no_rep$Time_body_out_sec_raw & summary_tracking_complete_no_rep$Time_body_out_sec_raw < 3000])
summary(fit.detect)

# Table with slope estimates of each linear fit
rbind(
  Censored = c(coef(fit.cen)[1], coef(fit.cen)[2]),
  OLS_Full = coef(fit.OLS),
  OLS_Detects = coef(fit.detect)
) %>% round(4)

# Plot the results
lineplot <- function() {
  
  # Add lines
  abline(coef(fit.cen)[1:2], col = 'Red', lwd = 2)
  abline(coef(fit.OLS), col = 'Blue', lty = 2, lwd = 2)
  abline(coef(fit.detect), col = rgb(0.2, 0.6, 0.2), lty = 3, lwd = 2)
  
  # Add legend
  legend(
    "topleft",
    legend = c(
      "Censored",
      "OLS",
      "OLS Detects"
    ),
    lwd = 3 * par("cex"),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 3),
    text.font = 1, bty = "n",
    pt.cex = 1, cex = 0.8, y.intersp = 1.3
  )
}

par(mfrow = c(1, 2))

plot(summary_tracking_complete_no_rep$Year_invasion, summary_tracking_complete_no_rep$Time_body_out_sec_raw, pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')
lineplot()

hist(summary_tracking_complete_no_rep$Time_body_out_sec_raw, breaks = 100, main = 'Censored Data')




# Visual exploration

plot(Visual_exp ~ Year_invasion, data = summary_tracking_complete_no_rep,  pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')

# Censored regression
fit.cen <- censReg(Visual_exp ~ Year_invasion, data = summary_tracking_complete_no_rep, 
                   left = 0, right = Inf)
summary(fit.cen)


# Get the standard error. It is pretty big (more than a thousand)
coef(fit.cen, logSigma = FALSE)


# Fit a line using OLS regression (linear regression)
fit.OLS <- lm(Visual_exp ~ Year_invasion, data = summary_tracking_complete_no_rep)
summary(fit.OLS)

# fit regression only of non-censored data
fit.detect <- lm(summary_tracking_complete_no_rep$Visual_exp[0 < summary_tracking_complete_no_rep$Visual_exp & summary_tracking_complete_no_rep$Visual_exp <= 3000] ~ summary_tracking_complete_no_rep$Year_invasion[0 < summary_tracking_complete_no_rep$Visual_exp & summary_tracking_complete_no_rep$Visual_exp <= 3000])
summary(fit.detect)

# Table with slope estimates of each linear fit
rbind(
  Censored = c(coef(fit.cen)[1], coef(fit.cen)[2]),
  OLS_Full = coef(fit.OLS),
  OLS_Detects = coef(fit.detect)
) %>% round(4)

# Plot the results
lineplot <- function() {
  
  # Add lines
  abline(coef(fit.cen)[1:2], col = 'Red', lwd = 2)
  abline(coef(fit.OLS), col = 'Blue', lty = 2, lwd = 2)
  abline(coef(fit.detect), col = rgb(0.2, 0.6, 0.2), lty = 3, lwd = 2)
  
  # Add legend
  legend(
    "topleft",
    legend = c(
      "Censored",
      "OLS",
      "OLS Detects"
    ),
    lwd = 3 * par("cex"),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 3),
    text.font = 1, bty = "n",
    pt.cex = 1, cex = 0.8, y.intersp = 1.3
  )
}

par(mfrow = c(1, 2))

plot(summary_tracking_complete_no_rep$Year_invasion, summary_tracking_complete_no_rep$Visual_exp, pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')
lineplot()

hist(summary_tracking_complete_no_rep$Visual_exp, breaks = 100, main = 'Censored Data')







plot(Visual_exp ~ Time_head_out_sec_raw, data = summary_tracking_complete_no_rep,  pch = 19, cex = 0.7, col = 'grey60', main = 'Censored Data')

fit.cen <- censReg(Visual_exp ~ Year_invasion + Temperature + Time_head_out_sec_raw, data = summary_tracking_complete_no_rep, 
                   left = 0, right = Inf)
summary(fit.cen)




hist(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), sd = sd(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))

# Head-out data follows log-normal distribution





# Example 2: Tobit model

install.packages("GGally")
install.packages("VGAM")
library(GGally)
library(VGAM)


ggpairs(summary_tracking_complete_no_rep[, c("Time_head_out_sec_raw", "Time_body_out_sec_raw", "Visual_exp")])

summary(m <- vglm(Time_head_out_sec_raw ~ Temperature + Time_body_out_sec_raw + Year_invasion + Inv_situation, tobit(Upper = 3000), data = summary_tracking_complete_no_rep))

summary(m <- vglm(Time_body_out_sec_raw ~ Year_invasion, tobit(Upper = 3000), data = summary_tracking_complete_no_rep))

summary(m <- vglm(Visual_exp ~ Year_invasion, tobit(Lower = 0), data = summary_tracking_complete_no_rep))

b <- coef(m)
se <- sqrt(diag(vcov(m)))

cbind(LL = b - qnorm(0.975) * se, UL = b + qnorm(0.975) * se)



summary_tracking_complete_no_rep$yhat <- fitted(m)[,1]
summary_tracking_complete_no_rep$rr <- resid(m, type = "response")
summary_tracking_complete_no_rep$rp <- resid(m, type = "pearson")[,1]

par(mfcol = c(2, 3))

with(summary_tracking_complete_no_rep, {
  plot(yhat, rr, main = "Fitted vs Residuals")
  qqnorm(rr)
  plot(yhat, rp, main = "Fitted vs Pearson Residuals")
  qqnorm(rp)
  plot(Time_head_out_sec_raw, rp, main = "Actual vs Pearson Residuals")
  plot(Time_head_out_sec_raw, yhat, main = "Actual vs Fitted")
})


(r <- with(summary_tracking_complete_no_rep, cor(yhat, Time_head_out_sec_raw)))

r^2


summary(m <- vglm(Visual_exp ~ Temperature + Time_head_out_sec_raw + Year_invasion + Inv_situation + Acc_time + Arena, tobit(Lower = 0), data = summary_tracking_complete_no_rep))



#####



##### Part 10: Data structure #####

# Head-out

par(mfrow = c(1,2))
hist((summary_tracking_complete_no_rep$Time_head_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep$Time_head_out_sec_raw)), sd = sd((summary_tracking_complete_no_rep$Time_head_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep$Time_head_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep$Time_head_out_sec_raw))

# Right censored data, delete censored data

hist((summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), sd = sd((summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))

# No normality. log transform data

hist(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), sd = sd(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test(log(summary_tracking_complete_no_rep_outside$Time_head_out_sec_raw))


# Body-out

hist((summary_tracking_complete_no_rep$Time_body_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep$Time_body_out_sec_raw)), sd = sd((summary_tracking_complete_no_rep$Time_body_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep$Time_body_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep$Time_body_out_sec_raw))

# Right censored data, delete censored data

hist((summary_tracking_complete_no_rep_outside$Time_body_out_sec_raw), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep_outside$Time_body_out_sec_raw)), sd = sd((summary_tracking_complete_no_rep_outside$Time_body_out_sec_raw))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep_outside$Time_body_out_sec_raw)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep_outside$Time_body_out_sec_raw))

# Normality


#Visual exploration

hist((summary_tracking_complete_no_rep$Visual_exp), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep$Visual_exp)), sd = sd((summary_tracking_complete_no_rep$Visual_exp))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep$Visual_exp)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep$Visual_exp))

# Left censored data, delete it

hist((summary_tracking_complete_no_rep_outside$Visual_exp), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep_outside$Visual_exp)), sd = sd((summary_tracking_complete_no_rep_outside$Visual_exp))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep_outside$Visual_exp)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep_outside$Visual_exp))

# No normality, log transoform

hist(log(summary_tracking_complete_no_rep_outside$Visual_exp), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean(log(summary_tracking_complete_no_rep_outside$Visual_exp)), sd = sd(log(summary_tracking_complete_no_rep_outside$Visual_exp))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density(log(summary_tracking_complete_no_rep_outside$Visual_exp)), main = "Gráfico de densidad de tus datos")

shapiro.test(log(summary_tracking_complete_no_rep_outside$Visual_exp))





##### Part 10.1: Residuals structure

par(mfrow = c(1,2))

hist((summary_tracking_complete_no_rep_outside$res_visual_svl_sex), prob = TRUE, main = "Histograma de tus datos")

curve(dnorm(x, mean = mean((summary_tracking_complete_no_rep_outside$res_visual_svl_sex)), sd = sd((summary_tracking_complete_no_rep_outside$res_visual_svl_sex))), 
      col = "blue", lwd = 2, add = TRUE)

plot(density((summary_tracking_complete_no_rep_outside$res_visual_svl_sex)), main = "Gráfico de densidad de tus datos")

shapiro.test((summary_tracking_complete_no_rep_outside$res_visual_svl_sex))

ggpairs(summary_tracking_complete_no_rep_outside[, c("res_head_cond", "res_body_cond", "res_visual_cond")])


pairs(~ res_head_svl + res_body_svl + res_visual_svl, 
      data = summary_tracking_complete_no_rep_outside,
      diag.panel = NULL)






#####







#####
