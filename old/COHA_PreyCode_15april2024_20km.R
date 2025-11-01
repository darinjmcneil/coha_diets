library(sp); library(raster); library(geosphere)
library(sf); library(ggplot2); library(tictoc)

# read in files
coha_all <- read.csv("./COHA_RawData.csv") # read in COHA raw data from working directory
prey_table <- read.csv("./PreyWeightTable_21jan2024.csv") # read in COHA prey weights from working directory
prey_table <- prey_table[ ,c(1,3,4,6)] # subset to get columns 1, 3, 4, and 6 (prey ID, min/max weight, and avian size class)

# merge coha_all with prey table to add all weights to identified prey
coha_all2 <- merge(x = coha_all, y = prey_table, by.x = "prey_desc", by.y = "prey_desc", all.x = TRUE)

# remove records that have no prey
coha_all2 <- subset(coha_all2, prey_desc != "no_prey")

# remove unidentified prey
coha_all2 <- subset(coha_all2, prey_desc != "?")
coha_all2 <- subset(coha_all2, prey_desc != "remove")

# remove unidentified avian
coha_all2 <- subset(coha_all2, prey_desc != "avian")

# remove unidentified mammalian
coha_all2 <- subset(coha_all2, prey_desc != "mammalian")

# add rownumber to coha_all2
coha_all2 <- cbind("row" = seq(1, nrow(coha_all2)), coha_all2)

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

# plotting squares
xwest <- us1@bbox[1,1] # use oh@bbox[1,1] for ohio runthrough; use us1@bbox[1,1] for USA runthrough
xeast <- us1@bbox[1,2]
ynorth <- us1@bbox[2,2]
ysouth <- us1@bbox[2,1]
width <- 10000
height <- 10000

# Note about prey size class replacements:
# On 15 April 2024, we tightened up the code a lot to add the
# prey replacements into the for() loops. This meant removing a lot
# of notes and such. See 1 April 2024 version that includes the old notes.

############################## Random weight assignment, random spatial subsampling, and replicating

# create blank list to hold resulting data.frames
list1 <- list()

for(i in 1:50){ # for each replicate (doing 50 runs through the data)
  #######################################################################
  # i = 34
  
  #######################################################################
  # assign random species IDs to small/medium/large avian prey items in random fashion
  
  ### "Small avian" replacement
  MysterySmalls <- subset(coha_all2, prey_desc == "small avian") # all unidentified small birds subsetted
  smallBirds <- subset(coha_all2, avian_size_class == "Small") # all identified small birds
  MysterySmallsUpd <- MysterySmalls # MysterySmallsUpd will contain updated species identities
  for(b in 1:nrow(MysterySmalls)){
    record_b <- MysterySmalls[b,] # get bth MysterySmall
    dist_df <- rbind(record_b, smallBirds) # we need to bind this mystery bird with
    coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
    crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
    d <- geosphere::distm(coords_sp)# calculate all pairwise distances
    smallBirds_b <- smallBirds # duplicate smalls for the purposes of the bth bird
    smallBirds_b$dist_to_record_b <- d[2:nrow(d), 1] # only need first column which is for record_b; also we don't need the first row b/c its the distance to itself
    smallBirds_b_10 <- subset(smallBirds_b, dist_to_record_b <= 10000)   # distance to nearest; all birds within 10km
    RandomCloseBird <- smallBirds_b_10[sample(1:nrow(smallBirds_b_10), 1),]  # parts of ifelse
    ClosestBird <- subset(smallBirds_b, dist_to_record_b == min(smallBirds_b$dist_to_record_b))
    ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird) # if closest bird is duplicated, select first record
    ifelse(nrow(smallBirds_b_10) > 0, NewRow <- RandomCloseBird, NewRow <- ClosestBird)# ifelse... if any birds are within 10km, sample a random one. If not, pick the closest
    NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
    NewRow[,3:8] <- record_b[,3:8] #
    NewRow$notes <- "semi-randomly assigned 'small' species"                           
    NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column
    NewRow$row <- record_b$row # retain its original row ID
    MysterySmallsUpd[b,] <- NewRow   # replace old row in MysterySmallsUpd with the new identity     
    print(paste0("small bird number ", b, " replaced with ", NewRow$prey_desc))  # print progress
  }
  print(paste0("replicate ", i, " small prey identities randomly assigned 🐦"))
  # there are 185 small birds
  
  ### "Medium avian" replacement
  MysteryMediums <- subset(coha_all2, prey_desc == "medium avian") # all unidentified medium birds subsetted
  mediumBirds <- subset(coha_all2, avian_size_class == "Medium") # all identified medium birds
  MysteryMediumsUpd <- MysteryMediums # MysteryMediumsUpd will contain updated species identities
  for(b in 1:nrow(MysteryMediums)){
    record_b <- MysteryMediums[b,] # get bth MysteryMedium
    dist_df <- rbind(record_b, mediumBirds) # we need to bind this mystery bird with
    coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
    crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
    d <- geosphere::distm(coords_sp)
    mediumBirds_b <- mediumBirds # duplicate mediums for the purposes of the bth bird
    mediumBirds_b$dist_to_record_b <- d[2:nrow(d), 1]
    mediumBirds_b_10 <- subset(mediumBirds_b, dist_to_record_b <= 10000) # all birds within 10km
    RandomCloseBird <- mediumBirds_b_10[sample(1:nrow(mediumBirds_b_10), 1),]
    ClosestBird <- subset(mediumBirds_b, dist_to_record_b == min(mediumBirds_b$dist_to_record_b))
    ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird)
    ifelse(nrow(mediumBirds_b_10) > 0, NewRow <- RandomCloseBird, NewRow <- ClosestBird)
    NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
    NewRow[,3:8] <- record_b[,3:8] #
    NewRow$notes <- "semi-randomly assigned 'medium' species"                           
    NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column                          
    NewRow$row <- record_b$row # retain its original row ID
    MysteryMediumsUpd[b,] <- NewRow
    print(paste0("medium bird number ", b, " replaced with ", NewRow$prey_desc))
  }
  print(paste0("replicate ", i, " medium prey identities randomly assigned 🐦🐦"))
  # there are 262 medium birds
  
  #### "Large avian" replacement
  MysteryLarges <- subset(coha_all2, prey_desc == "large avian") # all unidentified large birds subsetted
  largeBirds <- subset(coha_all2, avian_size_class == "Large") # all identified large birds
  MysteryLargesUpd <- MysteryLarges # MysteryLargesUpd will contain updated species identities
  for(b in 1:nrow(MysteryLarges)){
    record_b <- MysteryLarges[b,] # get bth MysteryLarge
    dist_df <- rbind(record_b, largeBirds) # we need to bind this mystery bird with
    coords_sp <- SpatialPoints(coords = data.frame("x" = dist_df$longitude, "y" = dist_df$latitude)) # create spatial object
    crs(coords_sp) <- CRS("+init=epsg:4269") # define CRS as NAD83
    d <- geosphere::distm(coords_sp)
    largeBirds_b <- largeBirds # duplicate larges for the purposes of the bth bird
    largeBirds_b$dist_to_record_b <- d[2:nrow(d), 1]
    largeBirds_b_10 <- subset(largeBirds_b, dist_to_record_b <= 10000) # all birds within 10km
    RandomCloseBird <- largeBirds_b_10[sample(1:nrow(largeBirds_b_10), 1),]
    ClosestBird <- subset(largeBirds_b, dist_to_record_b == min(largeBirds_b$dist_to_record_b))
    ifelse(nrow(ClosestBird)>1, ClosestBird <- ClosestBird[1,], ClosestBird <- ClosestBird)
    ifelse(nrow(largeBirds_b_10) > 0, NewRow <- RandomCloseBird, NewRow <- ClosestBird)
    NewRow$prey_desc <- paste0("r_", NewRow$prey_desc[1]) # add an "r_" to the species to indicate "randomly assigned"
    NewRow[,3:8] <- record_b[,3:8] #
    NewRow$notes <- "semi-randomly assigned 'large' species"                           
    NewRow <- NewRow[,1:ncol(NewRow)-1] # remove last column                          
    NewRow$row <- record_b$row # retain its original row ID
    MysteryLargesUpd[b,] <- NewRow
    print(paste0("large bird number ", b, " replaced with ", NewRow$prey_desc))
  }
  print(paste0("replicate ", i, " large prey identities randomly assigned 🐦🐦🐦"))
  # there are 147 large birds
  
  ### Recombine
  coha_all3 <- subset(coha_all2, prey_desc != "small avian" & # this subsets coha_all2 but REMOVES unidentified BIRDS
                        prey_desc != "medium avian" &
                        prey_desc != "large avian")
  coha_all3$min_weight <- as.numeric(coha_all3$min_weight) # making this numeric
  coha_all3$max_weight <- as.numeric(coha_all3$max_weight) # making this numeric
  coha_all3 <- subset(coha_all3, max_weight < 5400) # remove deer, raccoon, opossum, and cat this way
  coha_all3 <- rbind(coha_all3, MysterySmallsUpd, MysteryMediumsUpd, MysteryLargesUpd) # bind everything back together
  coha_all3$min_weight <- as.numeric(coha_all3$min_weight)
  coha_all3$max_weight <- as.numeric(coha_all3$max_weight)
  print(paste0("replicate ", i, " avian prey identities randomly assigned 🐦🐦🐦🐦")) 
  
  #######################################################################
  # assign weights to all prey items in random fashion
  coha_all3$weight <- 0  # blank column to hold randomly assigned prey weights
  for(j in 1:nrow(coha_all3)){
    # j = 1
    prey_j <- coha_all3[j,]
    low1 <- prey_j$min_weight # lower bound of prey distribution
    high1 <- prey_j$max_weight # upper bound of prey distribution
    mean1 <- mean(c(low1, high1)) # mean is the midpoint
    samp <- round(rnorm(n = 1, mean = mean1, sd = (high1-mean1)/2),0) # sample one point, round
    samp <- ifelse(samp < 0, mean1, samp)
    prey_j$weight <- samp # add new random weight to prey
    coha_all3[j,] <- prey_j # add column with new weight back to original data.frame
    #print(paste0("row ", j, " done. Prey ", prey_j$prey_desc, " assigned weight of ", samp, " g"))
  }
  print(paste0("replicate ", i, " prey weights randomly assigned ⚖️"))
  
  #######################################################################
  # create spatial object of all COHA observations
  Locs <- SpatialPoints(coords = data.frame("x" = coha_all3$longitude, "y" = coha_all3$latitude, "row" = coha_all3$row)) # create spatial object
  crs(Locs) <- CRS("+init=epsg:4269") # define CRS as NAD83
  Locs <- spTransform(Locs, crs(us1)) # reproject to match USA (or Ohio) raster
  
  # blank data frame to hold records
  RandomPreySample <- coha_all3[0,] # blank table that resembles coha_all3 without data
  
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
      #plot(poly, add = T, lwd = 5, border = "red")
      # extract points from the moving window polygon
      prey_table2 <- data.frame(Locs@coords)
      prey_table2 <- subset(prey_table2, x < max(poly[[1]][,1]) & x > min(poly[[1]][,1]))
      prey_table2 <- subset(prey_table2, y < max(poly[[1]][,2]) & y > min(poly[[1]][,2]))
      sampledRowNumber <- prey_table2[sample(1:nrow(prey_table2), size = 1),3]
      SelectedPrey <- subset(coha_all3, row == sampledRowNumber)
      RandomPreySample <- rbind(RandomPreySample, SelectedPrey)
    }

    percent_through <- (x - xwest) / (xeast - xwest) * 100
    percent_through <- round(percent_through, 0)
    print(paste0("Now ", percent_through, "% commplete through run ", i, "; ", nrow(RandomPreySample), " prey items so far."))
  }
  RandomPreySample[1:20,2]
  print(paste0("replicate ", i, " spatially balanced subset created 😎"))

  #######################################################################
  # add new data.frame to list1
  list1[[i]] <- RandomPreySample
}

plot(-1,-1, ylim = c(0, 0.005), xlim = c(0,500), xlab = "Prey Weight", ylab = "Probability",
     main = "Cooper's Hawk Prey Weights")

for(i in 1:length(list1)){
  RandomPreySample_i <- list1[[i]]
  c1 <- subset(RandomPreySample_i, weight < 500) # object "c1" is just for the histogram on the next line
  lines(density(c1$weight), lwd = 1, col = "black") # density plot line
}
rug(c1$weight, col = "red")

# run a single line with a rug at the bottom
plot(-1,-1, ylim = c(0, 0.005), xlim = c(0,500), xlab = "Prey Weight", ylab = "Probability")
RandomPreySample_i <- list1[[1]]
c1 <- subset(RandomPreySample_i, weight < 500) # object "c1" is just for the histogram on the next line
lines(density(c1$weight), lwd = 2, col = "black") # density plot line
rug(c1$weigh)

############################################ Proportions of prey

list2 <- list()

for(i in 1:length(list1)){
  # i = 1
  list_i <- list1[[i]] # isolate list i
  
  # remove the "random" assignments
  list_i <- list_i[!grepl("r_", list_i$prey_desc),] # remove unidentified prey items
  
  # remove the prey items that were not identified to species (they contain a slash "/")
  list_i <- list_i[!grepl("/", list_i$prey_desc),]
  
  # take a 50% subset
  list_i$use <- rbinom(n = nrow(list_i), size = 1, prob = 0.5) # assign random 50% subset
  list_i <- subset(list_i, use == 1)# take 50% subset
  
  # blank data.frame to hold list_i's proportions
  list_i_prey_props <- data.frame("species" = 0, "prop" = 0)
  
  # run through each prey item and add it to blank data.frame()
  for(s in 1:length(unique(list_i$prey_desc))){
    list_i_prey_s <- subset(list_i, prey_desc == unique(list_i$prey_desc)[s]) # isolate prey s
    prop_s <- 100 * nrow(list_i_prey_s)/nrow(list_i) # proportion = 100 x prey_s / all prey
    prop_s <- round(prop_s, 2) # round to 2 decimal places
    newrow <- c(unique(list_i$prey_desc)[s], prop_s)
    list_i_prey_props <- rbind(list_i_prey_props, newrow)
    list_i_prey_props$prop <- as.numeric(list_i_prey_props$prop)
  }
  
  list_i_prey_props <- list_i_prey_props[2:nrow(list_i_prey_props),] # delete blank first row
  list2[[i]] <- list_i_prey_props
}

############################################ Calculate means and 95% CIs for proportions

Forsplist <- coha_all2

# remove the "random" assignments and unidentified birds and slashes
Forsplist <- Forsplist[!grepl("/", Forsplist$prey_desc),] # remove unidentified prey items
Forsplist <- Forsplist[!grepl("avian", Forsplist$prey_desc),] # remove unidentified avian

# replace a few vague ones
Forsplist$prey_desc[grepl("norway rat", Forsplist$prey_desc)] <- "rat sp."
Forsplist$prey_desc[grepl("sylvilagus sp.", Forsplist$prey_desc)] <- "leporid sp." # sylvilagus becomes leporidae
Forsplist$prey_desc[grepl("eastern cottontail", Forsplist$prey_desc)] <- "leporid sp." # cottontail becomes leporidae
Forsplist$prey_desc[grepl("black-tailed jackrabbit", Forsplist$prey_desc)] <- "leporid sp." # cottontail becomes leporidae

splist <- unique(Forsplist$prey_desc)

PropData <- data.frame("species" = 0, "meanProp" = 0, "sdProp" = 0)

for(sp in 1:length(unique(Forsplist$prey_desc))){ # for each prey item...
  
  #sp = 9
  unique(Forsplist$prey_desc)[sp] # sp 53 is Eurasian Collared-dove
  
  sp_vals <- c()
  
  for(l in 1:length(list2)){ # obtain all of the proportions from list2
    prop_l <- subset(list2[[l]], species == unique(Forsplist$prey_desc)[sp]) 
    prop_l <- as.numeric(prop_l[2])
    sp_vals <- c(sp_vals, prop_l)
  }
  
  # then calculate the mean and SD for each
  speciesName = unique(Forsplist$prey_desc)[sp]
  meanProportion = round(mean(sp_vals, na.rm = T),2)
  sdProportion = round(sd(sp_vals, na.rm = T),2)
  
  # add them to data.frame()
  PropData <- rbind(PropData, c(speciesName, meanProportion, sdProportion))
}

PropData <- PropData[2:nrow(PropData),] # remove blank row
PropData <- PropData[complete.cases(PropData),] # remove NAs
PropData$meanProp <- as.numeric(PropData$meanProp) # make means numeric
PropData$sdProp <- as.numeric(PropData$sdProp) # make sds numeric
PropData <- PropData[order(PropData$meanProp, decreasing = TRUE),]
PropData$seProp <- PropData$sdProp/sqrt(length(list2))
PropData$ciProp <- 1.96 * PropData$seProp
PropData[,4:5] <- round(PropData[,4:5], 2)
PropData


