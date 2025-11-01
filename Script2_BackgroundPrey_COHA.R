library(rgbif); library(sf); library(sp); library(terra);library(eSDM)
library(dplyr); library(utils); library(dwctaxon); library(readr); library(tibble)
library(tidyr); library(maps); library(stringr)

#bringing in iNat background data (one class at a time)
allaves <- data.frame(read.csv("C:/Users/User/Desktop/COHA_code/all_bg_prey/allaves.csv"))
allmamm <- data.frame(read.csv("C:/Users/User/Desktop/COHA_code/all_bg_prey/allmamm.csv"))
allrept <- data.frame(read.csv("C:/Users/User/Desktop/COHA_code/all_bg_prey/allrept.csv"))

#remove NAs in iNat background data
allaves <- tidyr::drop_na(allaves)
allrept <- tidyr::drop_na(allrept)
allmamm <- tidyr::drop_na(allmamm)

#US map
us <- terra::vect("C:/Users/User/Desktop/COHA_code/states/States_shapefile.shp")
us <- us[us$State_Code != 'AK' & us$State_Code != 'HI',] # remove AK and HI
us1 <- terra::project(us, crs("+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")) # reproject in lambert conic projection; makes creating grids (in km) easier
plot(us1, axes = T); plot(us[us$State_Code == 'AL', ])

#generate spatial objects for iNat background data
# birds
birdpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
birdpt <- terra::vect(birdpt)
terra::crs(birdpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
birdpt <- terra::project(birdpt, crs(us1))
# reptiles
reptpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
reptpt = terra::vect(reptpt)
terra::crs(reptpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
reptpt = terra::project(reptpt, crs(us1))
# mammals
mammpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
mammpt = terra::vect(mammpt)
terra::crs(mammpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
mammpt = terra::project(mammpt, crs(us1))

#####################

#read in all COHA depredation observations (THINNED records of hawks eating prey)
list1 <- readRDS("./list1.rds")
coha <- list1[[1]]

#make coha coordinates named x and y for later ease
names(coha)[names(coha) == "latitude"] = "y"
names(coha)[names(coha) == "longitude"] = "x"

#remove s, m, l avian (these have "r_", at the begining from Script1)
coha <- coha[!grepl("r_", coha$prey_desc), ]

#subset prey classes
cohaaves <- subset(coha, prey_class = "aves")
cohamammals <- subset(coha, prey_class = "mammalia")
cohaherp <- subset(coha, prey_class = "sauria")

#generate planar projection for thinned coha depredation points
#cohacoords <- dplyr::select(coha, x, y) #subset coordinates
#Locs <- SpatialPoints(coords = data.frame("x" = cohacoords$x, "y" = cohacoords$y)) # create spatial object
#Locs <- terra::vect(Locs)
#terra::crs(Locs) <- CRS("+init=epsg:4269") # define CRS as NAD83
#Locs <- terra::project(Locs, crs(us1)) # reproject to match map of the US

### better than above chunk?
library(terra)
Locs <- vect(coha, geom = c("x", "y"), crs = "EPSG:4269")
Locs <- project(Locs, crs(us1))

completeprey_df = data.frame("coha_record" = 0,
                             "coha_date" = 0,
                             "state" = 0,
                             "coha_prey" = 0,
                             "prey_class" = 0,
                             "rand_prey" = 0,
                             "rand_prey_record" = 0,
                             "rand_date" = 0,
                             "state" = 0,
                             "NA" = 0)

for(i in 1:length(us$State_Code)){ #for every state
  #i = 1
  
  # isolate the name of state i
  state_i <- us$State_Code[i]
  
  # obtain "fat" (buffered) state i boundary and plot it
  str3 <- paste0("st_i_buff <- terra::vect('C:/Users/User/Desktop/COHA_code/StateStuff/statebuff/", state_i, "_buffered.shp')")
  eval(parse(text = str3)) #actually runs str0 as code
  # plot(st_i_buff, add = T, lwd = 3) # plotting turned off for speed purposes
  
  # crop actual predation records to the "fat" (buffered) state i boundary and plot it
  st_i_coha <- terra::crop(Locs, us1[i])
  # plot(st_i_coha, add = T, col = "red")  # plotting turned off for speed purposes
  
  if (length(st_i_coha) == 0) {
    next   # skips to the next iteration in a for loop
  } else {
    
    # obtain the background prey from the "stateiNatbirds" folder for state i
    str2 <- paste0("stateiNatbirds <- terra::vect('C:/Users/User/Desktop/COHA_code/StateStuff/stateiNatbirds/", state_i, "_iNatbirds.shp')")
    eval(parse(text = str2)) #actually runs str0 as code
    plot(stateiNatbirds, main = paste0(state_i)) # plotting turned off for speed purposes
    
    # This for() loop runs through all the COHA depredation records for state i 
    # and obtains 1000 prey items for each event  
    
    for(j in 1:length(st_i_coha)){ #loops through each coha "j" within state "i"
      #j = 1
      #tic()
      buff1 = terra::buffer(st_i_coha[j], 25000)
      plot(buff1, add = T, border = "lightcoral", lwd = 2)
      cohaj_potprey = terra::crop(stateiNatbirds, buff1) #potential prey for coha j
      
      if (length(length(cohaj_potprey)) == 0) {
        next   # skips to the next iteration in a for loop
      } else {plot(cohaj_potprey, add = T, col = "cyan")}
      
      preysamp_j = sample(cohaj_potprey, 1000, replace = T) #sample of pot prey
      plot(preysamp_j, add = T, col = "purple")
      cohaj_randpreyitems <- data.frame("coha_record" = st_i_coha@ptr$df$values()$url[j],
                                        "coha_date" = st_i_coha@ptr$df$values()$observed_on[j],
                                        "state" = state_i,
                                        "coha_prey" = st_i_coha@ptr$df$values()$prey_desc[j],
                                        "prey_class" = "Aves",
                                        "rand_prey" = preysamp_j@ptr$df$values()$species,
                                        "rand_prey_record" = preysamp_j@ptr$df$values()$url,
                                        "rand_date" = preysamp_j@ptr$df$values()$date,
                                        "state.1" = preysamp_j@ptr$df$values()$state,
                                        "NA." = "NA")
      names(cohaj_randpreyitems) = c("coha_record", "coha_date", "state", "coha_prey",
                                     "prey_class", "rand_prey", "rand_prey_record",
                                     "rand_date", "state.1", "NA.")
      completeprey_df = rbind(completeprey_df, cohaj_randpreyitems)
      print(paste0("state number ", i, " (", state_i, ") ", "bird ", j, " of ", length(st_i_coha), " done"))
    
      } # ends for() loop that runs through all records within a state
    } # ends if() statement that checks for records within a state
  }# ends for() loop that runs through each state

#
#
#
#
#
#
#
#
#
#
#
#
#
#


###########################
#random selection
#birds
plot(us1)
birdselect <- data.frame(row.names = "species", "latitude", "longitude")

for(i in 1:5){
  buff1 = terra::buffer(Locs[i], 10000)
  wkt = as.character(buff1)
  plot(buff1, add = T)
  birdclip = terra::crop(birdpt, buff1)
  plot(birdclip, add = T)
  randbird = sample(birdclip, 1)
  plot(randbird, add = T, col = "red")
}

rbind(randbird, birdselect)

#reptiles
plot(us1)
reptselect = data.frame(row.names = "species", "latitude", "longitude")
for(i in 1:5){
  buff1 = terra::buffer(Locs[i], 50000)
  wkt = as.character(buff1)
  plot(buff1, add = T)
  reptclip = terra::crop(reptpt, buff1)
  plot(reptclip, add = T)
  randrept = sample(reptclip, 1)
  plot(randrept, add = T, col = "red")
}

#mammals
plot(us1)
mammselect = data.frame(row.names = "species", "latitude", "longitude")
for(i in 1:5){
  buff1 = terra::buffer(Locs[i], 50000)
  wkt = as.character(buff1)
  plot(buff1, add = T)
  mammclip = terra::crop(mammpt, buff1)
  plot(mammclip, add = T)
  randmamm = sample(mammclip, 1)
  plot(randmamm, add = T, col = "red")
}
###################################
#download full US dataset
#trim all the extra columns and everything- just have spp inat record and latlong
#turn them into a spatial object and set it aside
#take COHA record 100 and buffer by 25km then clip bird dataframe by just COHA buffer


#creating a polygon around prey points

#birdpoly = pts2poly_centroids(Locs, 100)
#birdpoly = st_as_sf(birdpoly$geometry)
#birdpoly = as_Spatial(birdpoly)
#mammalpoly = pts2poly_centroids(cohamammals, 25, crs = crs(us1))
#herppoly = pts2poly_centroids(cohaherp, 25, crs = crs(us1))

#making that polygon into WKT

###############
#creating keys
#dataset key used beow ("50c9509d-22c7-4a22-a47d-8c48425ef4a7") is key for iNat data
#birdkey = name_backbone(name = "Aves"); birdkey = birdkey$usageKey #key for birds
#mammalkey = name_backbone(name = "Mammalia"); mammalkey = mammalkey$usageKey #key for mammals
#reptilekey = name_backbone(name = "Squamata"); reptilekey = reptilekey$usageKey #key for reptiles

##############
#occupancy search for bird record

#ex1 = occ_search(taxonKey = birdkey,
#                 datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7",
#                 hasCoordinate = T, geometry = birdpoly[[3]][[1]][[1]],
#)

#ex1 = ex1$data
#ex1= ex1$scientificName
#View(ex1)

#samp = spsample(birdpoly[1],n=1,"random")
#samp = st_as_sf(samp)
#st_coordinates(samp)

#plot point, nearest record to that point

#for(i in 1:50) {
#  plot(birdpoly[i])
#  samp = spsample(birdpoly[i], n=1, "random")
#  samp = st_as_sf(samp)
#  samp = st_coordinates(samp)
#  samp = bbox(samp)
#  samp = as.data.frame(samp)
#  samp = gbif_bbox2wkt(samp$min[1], samp$min[2], samp$max[1], samp$max[2])
#}

#occ_search(taxonKey = birdkey,
#           datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7",
#           hasCoordinate = T, geometry = cplanar)
