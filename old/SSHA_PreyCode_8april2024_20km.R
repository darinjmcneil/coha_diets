library(sp); library(raster); library(rgeos); library(geosphere)
library(sf); library(ggplot2); library(tictoc)

# read in files
ssha_all <- read.csv("./SSHA_RawData.csv") # read in SSHA raw data from working directory
prey_table <- read.csv("./SSHA_PreyWeightTable_8april2024.csv") # read in COHA prey weights from working directory
prey_table <- prey_table[ ,c(1,3,4,6)] # subset to get columns 1, 3, 4, and 6 (prey ID, min/max weight, and avian size class)

# merge coha_all with prey table to add all weights to identified prey
ssha_all2 <- merge(x = ssha_all, y = prey_table, by.x = "prey_desc", by.y = "item", all.x = TRUE)

# remove records that have no prey
ssha_all2 <- subset(ssha_all2, prey_desc != "")

# remove unidentified prey
ssha_all2 <- subset(ssha_all2, prey_desc != "?")
ssha_all2 <- subset(ssha_all2, prey_desc != "remove")

# remove unidentified avian
ssha_all2 <- subset(ssha_all2, prey_desc != "avian")

# remove records that lack coordinates
ssha_all2 <- ssha_all2[!(is.na(ssha_all2$latitude)), ]
ssha_all2 <- subset(ssha_all2, latitude > 25) # remove tropical records
ssha_all2 <- subset(ssha_all2, latitude < 50) # remove arctic records

# add rownumber to ssha_all2
ssha_all2 <- cbind("row" = seq(1, nrow(ssha_all2)), ssha_all2)

##########################################################
########################################################## "Small avian" replacement
##########################################################

### Assign random species to each "small prey"
MysterySmalls <- subset(ssha_all2, prey_desc == "small avian") # all unidentified small birds subsetted
smallBirds <- subset(ssha_all2, avian.size.class == "Small") # all identified small birds

# https://hub.arcgis.com/datasets/1b02c87f62d24508970dc1a6df80c98e/explore?location=38.753686%2C-99.481986%2C3.89
# https://stackoverflow.com/questions/21977720/r-finding-closest-neighboring-point-and-number-of-neighbors-within-a-given-rad
#us <- raster::shapefile("C:/Users/User/Desktop/COHA_code/States_shapefile.shp") # Read in shapefile of US
#us <- us[us$State_Code != 'AK' & us$State_Code != 'HI',] # remove AK and HI
#us1 <- spTransform(us, CRS("+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"))# reproject in lambert conic projection; makes creating grids (in km) easier
#plot(us1)

################################################################################

MysterySmallsUpd <- MysterySmalls # MysterySmallsUpd will contain updated species identities

for(b in 1:nrow(MysterySmalls)){
  # b = 177
  record_b <- MysterySmalls[b,] # get bth MysterySmall
  dist_df <- rbind(record_b, smallBirds) # we need to bind this mystery bird with
  # the identified birds for the little spatial analysis below
  
  # convert raw coordinates to spatial data
  coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
  crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
  
  # calculate all pairwise distances
  d <- geosphere::distm(coords_sp)
  
  # only need first column which is for record_b
  # also we don't need the first row b/c its the distance to itself
  smallBirds_b <- smallBirds # duplicate smalls for the purposes of the bth bird
  smallBirds_b$dist_to_record_b <- d[2:nrow(d), 1]
  
  # distance to nearest 
  smallBirds_b_10 <- subset(smallBirds_b, dist_to_record_b <= 10000) # all birds within 10km
  #print(paste0("bird number ", b, " has ", nrow(smallBirds_b_10), " obs nearby"))
  
  # parts of ifelse
  RandomCloseBird <- smallBirds_b_10[sample(1:nrow(smallBirds_b_10), 1),]
  ClosestBird <- subset(smallBirds_b, dist_to_record_b == min(smallBirds_b$dist_to_record_b))
  
  # if closest bird is duplicated, select first record
  ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird)
  
  # ifelse... if any birds are within 10km, sample a random one. If not, pick the closest
  ifelse(nrow(smallBirds_b_10) > 0, 
         NewRow <- RandomCloseBird, 
         NewRow <- ClosestBird)
  NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
  NewRow[,3:8] <- record_b[,3:8] #
  NewRow$notes <- "semi-randomly assigned 'small' species"                           
  NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column
  NewRow$row <- record_b$row # retain its original row ID
  
  # replace old row in MysterySmallsUpd with the new identity                        
  MysterySmallsUpd[b,] <- NewRow
  
  # print progress
  print(paste0("small bird number ", b, " replaced with ", NewRow$prey_desc))
}

##########################################################
########################################################## "Medium avian" replacement
##########################################################

### Assign random species to each "medium prey"
MysteryMediums <- subset(ssha_all2, prey_desc == "medium avian") # all unidentified medium birds subsetted
mediumBirds <- subset(ssha_all2, avian.size.class == "Medium") # all identified medium birds

################################################################################

MysteryMediumsUpd <- MysteryMediums # MysteryMediumsUpd will contain updated species identities

for(b in 1:nrow(MysteryMediums)){
  # b = 177
  record_b <- MysteryMediums[b,] # get bth MysteryMedium
  dist_df <- rbind(record_b, mediumBirds) # we need to bind this mystery bird with
  # the identified birds for the little spatial analysis below
  
  # convert raw coordinates to spatial data
  coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
  crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
  
  # calculate all pairwise distances
  d <- geosphere::distm(coords_sp)
  
  # only need first column which is for record_b
  # also we don't need the first row b/c its the distance to itself
  mediumBirds_b <- mediumBirds # duplicate mediums for the purposes of the bth bird
  mediumBirds_b$dist_to_record_b <- d[2:nrow(d), 1]
  
  # distance to nearest 
  mediumBirds_b_10 <- subset(mediumBirds_b, dist_to_record_b <= 10000) # all birds within 10km
  #print(paste0("bird number ", b, " has ", nrow(mediumBirds_b_10), " obs nearby"))
  
  # parts of ifelse
  RandomCloseBird <- mediumBirds_b_10[sample(1:nrow(mediumBirds_b_10), 1),]
  ClosestBird <- subset(mediumBirds_b, dist_to_record_b == min(mediumBirds_b$dist_to_record_b))
  
  # if closest bird is duplicated, select first record
  ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird)
  
  # ifelse... if any birds are within 10km, sample a random one. If not, pick the closest
  ifelse(nrow(mediumBirds_b_10) > 0, 
         NewRow <- RandomCloseBird, 
         NewRow <- ClosestBird)
  NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
  NewRow[,3:8] <- record_b[,3:8] #
  NewRow$notes <- "semi-randomly assigned 'medium' species"                           
  NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column                          
  NewRow$row <- record_b$row # retain its original row ID
  
  # replace old row in MysteryMediumsUpd with the new identity                        
  MysteryMediumsUpd[b,] <- NewRow
  
  # print progress
  print(paste0("medium bird number ", b, " replaced with ", NewRow$prey_desc))
}

##########################################################
########################################################## "Large avian" replacement
##########################################################

### Assign random species to each "large prey"
MysteryLarges <- subset(ssha_all2, prey_desc == "large avian") # all unidentified large birds subsetted
largeBirds <- subset(ssha_all2, avian.size.class == "Large") # all identified large birds

################################################################################

MysteryLargesUpd <- MysteryLarges # MysteryLargesUpd will contain updated species identities

for(b in 1:nrow(MysteryLarges)){
  # b = 177
  record_b <- MysteryLarges[b,] # get bth MysteryLarge
  dist_df <- rbind(record_b, largeBirds) # we need to bind this mystery bird with
  # the identified birds for the little spatial analysis below
  
  # convert raw coordinates to spatial data
  coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
  crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
  
  # calculate all pairwise distances
  d <- geosphere::distm(coords_sp)
  
  # only need first column which is for record_b
  # also we don't need the first row b/c its the distance to itself
  largeBirds_b <- largeBirds # duplicate larges for the purposes of the bth bird
  largeBirds_b$dist_to_record_b <- d[2:nrow(d), 1]
  
  # distance to nearest 
  largeBirds_b_10 <- subset(largeBirds_b, dist_to_record_b <= 10000) # all birds within 10km
  #print(paste0("bird number ", b, " has ", nrow(largeBirds_b_10), " obs nearby"))
  
  # parts of ifelse
  RandomCloseBird <- largeBirds_b_10[sample(1:nrow(largeBirds_b_10), 1),]
  ClosestBird <- subset(largeBirds_b, dist_to_record_b == min(largeBirds_b$dist_to_record_b))
  
  # if closest bird is duplicated, select first record
  ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird)
  
  # ifelse... if any birds are within 10km, sample a random one. If not, pick the closest
  ifelse(nrow(largeBirds_b_10) > 0, 
         NewRow <- RandomCloseBird, 
         NewRow <- ClosestBird)
  NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
  NewRow[,3:8] <- record_b[,3:8] #
  NewRow$notes <- "semi-randomly assigned 'large' species"                           
  NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column                          
  NewRow$row <- record_b$row # retain its original row ID
  
  # replace old row in MysteryLargesUpd with the new identity                        
  MysteryLargesUpd[b,] <- NewRow
  
  # print progress
  print(paste0("large bird number ", b, " replaced with ", NewRow$prey_desc))
}

##########################################################
########################################################## Recombining data 
##########################################################

# pull the randoms in with their newly assigned pseudo-species to create one
# complete dataset where every observation has a prey (or pseudo-prey) for which
# weight can be calculated
# NOTE about mammals: I don't think we can have pseudo-mammals b/c we don't have
# mammal size classes (e.g., small mammal, medium mammal)
# So, in other words, mammalian prey stays IN this analysis, but ONLY mammals
# that we identified to species; a random, unidentified mammal is not used.

ssha_all3 <- subset(ssha_all2, prey_desc != "small avian" & # this subsets ssha_all2 but REMOVES unidentified BIRDS
                               prey_desc != "medium avian" &
                               prey_desc != "large avian")
ssha_all3$min.weight..g. <- as.numeric(ssha_all3$min.weight..g.) # making this numeric
ssha_all3$max.weight..g. <- as.numeric(ssha_all3$max.weight..g.) # making this numeric
ssha_all3 <- subset(ssha_all3, max.weight..g. < 5400) # remove deer, raccoon, opossum, and cat this way
hist(ssha_all3$min.weight..g., breaks = 100)

# bind everything back together
ssha_all3 <- rbind(ssha_all3, MysterySmallsUpd, MysteryMediumsUpd, MysteryLargesUpd)
ssha_all3$min_weight <- as.numeric(ssha_all3$min.weight..g.) # making this numeric
ssha_all3$max_weight <- as.numeric(ssha_all3$max.weight..g.) # making this numeric

# remove extra weight columns
ssha_all3 <- dplyr::select(ssha_all3, -min.weight..g., -max.weight..g.)

############################## Read in shapefiles and create spatial objects

# download shapefile from ESRI
#https://hub.arcgis.com/datasets/1b02c87f62d24508970dc1a6df80c98e/explore?location=38.753686%2C-99.481986%2C3.89
us <- raster::shapefile("./States_shapefile.shp")
us <- us[us$State_Code != 'AK' & us$State_Code != 'HI',] # remove AK and HI
plot(us)
crs(us)

# reproject in lambert conic projection; makes creating grids (in km) easier
lam <- "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"
us1 <- spTransform(us, CRS(lam))
plot(us1, axes = T)

# testing with Ohio since it has a relatively even spread of points
#oh <- us1[us1$State_Code == 'OH',] # subset to just Ohio
#crs(oh) # identify CRS
#plot(oh, axes = T) # plot

# plotting squares
xwest <- us1@bbox[1,1] # use oh@bbox[1,1] for ohio runthrough; use us1@bbox[1,1] for USA runthrough
xeast <- us1@bbox[1,2]
ynorth <- us1@bbox[2,2]
ysouth <- us1@bbox[2,1]
width <- 10000
height <- 10000

############################## Random weight assignment, random spatial subsampling, and replicating

# create blank list to hold resulting data.frames
list1 <- list()

for(i in 1:50){
  #######################################################################
  # i = 33
  
  #######################################################################
  # assign weights to all prey items in random fashion
  ssha_all3$weight <- 0  # blank column to hold randomly assigned prey weights
  for(j in 1:nrow(ssha_all3)){
    # j = 1
    prey_j <- ssha_all3[j,]
    low1 <- prey_j$min_weight # lower bound of prey distribution
    high1 <- prey_j$max_weight # upper bound of prey distribution
    mean1 <- mean(c(low1, high1)) # mean is the midpoint
    samp <- round(rnorm(n = 1, mean = mean1, sd = (high1-mean1)/2),0) # sample one point, round
    samp <- ifelse(samp < 0, mean1, samp)
    prey_j$weight <- samp # add new random weight to prey
    ssha_all3[j,] <- prey_j # add column with new weight back to original data.frame
    #print(paste0("row ", j, " done. Prey ", prey_j$prey_desc, " assigned weight of ", samp, " g"))
  }
  print(paste0("replicate ", i, " prey weights randomly assigned 👍"))
  
  #######################################################################
  # create spatial object of all ssha observations
  Locs <- SpatialPoints(coords = data.frame("x" = ssha_all3$longitude, "y" = ssha_all3$latitude, "row" = ssha_all3$row)) # create spatial object
  crs(Locs) <- CRS("+init=epsg:4269") # define CRS as NAD83
  Locs <- spTransform(Locs, crs(us1)) # reproject to match USA (or Ohio) raster
  
  # blank data frame to hold records
  RandomPreySample <- ssha_all3[0,] # blank table that resembles ssha_all3 without data
  
  # for() loop that cycles through all Cooper's Hawk observations and extracts random ones
  for(x in seq(xwest, xeast, by = 20000)){ # for each value of xwest, by 10km chunks
    # x = 591507.1
    for(y in seq(ysouth, ynorth, by = 20000)){ # for each value of ynorth, by 10km chunks
      # y = 4743980
      poly <- sf::st_polygon(
        list(
          cbind(
            c(x - width, x + width, x + width, x - width, x - width),
            c(y + height, y + height, y - height, y - height, y + height))
        )
      )
      # extract points from the moving window polygon
      prey_table2 <- data.frame(Locs@coords)
      prey_table2 <- subset(prey_table2, x < max(poly[[1]][,1]) & x > min(poly[[1]][,1]))
      prey_table2 <- subset(prey_table2, y < max(poly[[1]][,2]) & y > min(poly[[1]][,2]))
      sampledRowNumber <- prey_table2[sample(1:nrow(prey_table2), size = 1),3]
      SelectedPrey <- subset(ssha_all3, row == sampledRowNumber)
      RandomPreySample <- rbind(RandomPreySample, SelectedPrey)
    }
  }
  print(paste0("replicate ", i, " spatially balanced subset created 😎"))

  #######################################################################
  # add new data.frame to list1
  list1[[i]] <- RandomPreySample
}

plot(-1,-1, ylim = c(0, 0.015), xlim = c(0,500), xlab = "Prey Weight", ylab = "Probability",
     main = "Sharp-shinned Hawk Prey Weights")

for(i in 1:length(list1)){
  RandomPreySample_i <- list1[[i]]
  c1 <- subset(RandomPreySample_i, weight < 500) # object "c1" is just for the histogram on the next line
  lines(density(c1$weight), lwd = 1, col = "black") # density plot line
}
rug(c1$weigh, col = "red")


# run a single line with a rug at the bottom
plot(-1,-1, ylim = c(0, 0.005), xlim = c(0,500), xlab = "Prey Weight", ylab = "Probability")
RandomPreySample_i <- list1[[1]]
c1 <- subset(RandomPreySample_i, weight < 500) # object "c1" is just for the histogram on the next line
lines(density(c1$weight), lwd = 2, col = "black") # density plot line
rug(c1$weigh)
