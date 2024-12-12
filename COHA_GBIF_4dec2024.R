library(rgbif)
library(sf)
library(sp)
library(terra)
library(eSDM)
library(dplyr)
library(utils)
library(dwctaxon)
library(readr)
library(tibble)
library(tidyr)

#cohacoords = subset(coha, select = latitude:longitude) #isolating lat and long

#################
##DwCA DOIs

#aves: 10.15468/dl.ucn73z
#reptiles: 10.15468/dl.j2k7yz
#mammals: 10.15468/dl.ewaavu

##################
##bring in data

#birds
#occ_download(pred("taxonKey", 212), pred("country", "US"), 
#             pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
#             format = "DWCA", user = "wgibsonky", pwd = "Wjg742611!", 
#             email = "wgibsonky@gmail.com")

#mammals
#occ_download(pred("taxonKey", 359), pred("country", "US"), 
#             pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
#             format = "DWCA", user = "wgibsonky", pwd = "Wjg742611!", 
#             email = "wgibsonky@gmail.com")

#reptiles (squamata)
#occ_download(pred("taxonKey", 11592253), pred("country", "US"), 
#             pred("datasetKey", "50c9509d-22c7-4a22-a47d-8c48425ef4a7"),
#             format = "DWCA", user = "wgibsonky", pwd = "Wjg742611!", 
#             email = "wgibsonky@gmail.com")


#read in all 1628 prey observations
coha = read.csv("./cohaprey.csv")

#bringing in iNat data

#data = read_tsv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/dwcabird/occurrence.txt")
#data = read_tsv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/dwcarept/occurrence.txt")
#data = read_tsv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/dwcamamm/occurrence.txt")

#data = dplyr::select(data, species, decimalLatitude, decimalLongitude)

allaves = data.frame(read.csv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/allaves.csv"))
allrept = data.frame(read.csv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/allrept.csv"))
allmamm = data.frame(read.csv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/coha_diets/allmamm.csv"))

###################
##data processing

#remove unneeded rows in coha
coha = dplyr::select(coha, -X, -observed_on, -url, -prey_Y.N, -notes, -observer,
                     -min_weight, -max_weight, -avian_size_class)

#make coha coordinates named x and y for later ease
names(coha)[names(coha) == "latitude"] = "y"
names(coha)[names(coha) == "longitude"] = "x"

#remove s, m, l avian
coha = subset(coha, prey_desc != "small avian")
coha = subset(coha, prey_desc != "medium avian")
coha = subset(coha, prey_desc != "large avian")

#subset prey classes
cohaaves = subset(coha, prey_class = "aves")
cohamammals = subset(coha, prey_class = "mammalia")
cohaherp = subset(coha, prey_class = "sauria")

#remove nas in inat data
allaves = tidyr::drop_na(allaves)
allrept = tidyr::drop_na(allrept)
allmamm = tidyr::drop_na(allmamm)

###################
#US map
#https://hub.arcgis.com/datasets/1b02c87f62d24508970dc1a6df80c98e/explore?location=38.753686%2C-99.481986%2C3.89
us = terra::vect("C:/Users/wgibs/OneDrive/Desktop/coha_SSHA Diet/States_shapefile.shp")
us = us[us$State_Code != 'AK' & us$State_Code != 'HI',] # remove AK and HI
plot(us)
crs(us)

# reproject in lambert conic projection; makes creating grids (in km) easier
lam = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"
us1 = terra::project(us, crs(lam))
plot(us1, axes = T)

################
#generate planar projection for points

cohacoords = dplyr::select(coha, x, y) #subset coordinates

Locs <- SpatialPoints(coords = data.frame("x" = cohacoords$x, "y" = cohacoords$y)) # create spatial object
Locs = terra::vect(Locs)
terra::crs(Locs) <- CRS("+init=epsg:4269") # define CRS as NAD83
Locs = terra::project(Locs, crs(us1))

#generate spatial objects for iNat data

birdpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
birdpt = terra::vect(birdpt)
terra::crs(birdpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
birdpt = terra::project(birdpt, crs(us1))

reptpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
reptpt = terra::vect(reptpt)
terra::crs(reptpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
reptpt = terra::project(reptpt, crs(us1))

mammpt <- SpatialPoints(coords = data.frame("x" = allaves$decimalLongitude, "y" = allaves$decimalLatitude)) # create spatial object
mammpt = terra::vect(mammpt)
terra::crs(mammpt) <- CRS("+init=epsg:4269") # define CRS as NAD83
mammpt = terra::project(mammpt, crs(us1))


#plotting 1 point

#for(i in 1:nrow(Locs)){
#  plot(us1)
#  plot(Locs[i])
#  buff1 = terra::buffer(Locs[i], 100000)
#}


plot(us1)
plot(Locs[100], add = T)
buff1 = terra::buffer(Locs[100], 100000)
wkt = as.character(buff1)
class(wkt)
plot(buff1, add = T)

as.data.frame(buff1)


#download full US dataset
#trim all the extra columns and everything- just have spp inat record and latlong
#turn them into a spatial object and set it aside
#take cCOHA record 100 and buffer by 25km then clip bird dataframe by just COHA buffer


#creating a polygon around prey points

birdpoly = pts2poly_centroids(Locs, 100)
birdpoly = st_as_sf(birdpoly$geometry)
birdpoly = as_Spatial(birdpoly)
mammalpoly = pts2poly_centroids(cohamammals, 25, crs = crs(us1))
herppoly = pts2poly_centroids(cohaherp, 25, crs = crs(us1))

#making that polygon into WKT

###############
#creating keys
#dataset key used beow ("50c9509d-22c7-4a22-a47d-8c48425ef4a7") is key for iNat data
birdkey = name_backbone(name = "Aves"); birdkey = birdkey$usageKey #key for birds
mammalkey = name_backbone(name = "Mammalia"); mammalkey = mammalkey$usageKey #key for mammals
reptilekey = name_backbone(name = "Squamata"); reptilekey = reptilekey$usageKey #key for reptiles

##############
#occupancy search for bird record

#ex1 = occ_search(taxonKey = birdkey,
#                 datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7",
#                 hasCoordinate = T, geometry = birdpoly[[3]][[1]][[1]],
#)

#ex1 = ex1$data
#ex1= ex1$scientificName
#View(ex1)

samp = spsample(birdpoly[1],n=1,"random")
samp = st_as_sf(samp)
st_coordinates(samp)

#plot point, nearest record to that point

for(i in 1:50) {
  plot(birdpoly[i])
  samp = spsample(birdpoly[i], n=1, "random")
  samp = st_as_sf(samp)
  samp = st_coordinates(samp)
  samp = bbox(samp)
  samp = as.data.frame(samp)
  samp = gbif_bbox2wkt(samp$min[1], samp$min[2], samp$max[1], samp$max[2])
}

occ_search(taxonKey = birdkey,
           datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7",
           hasCoordinate = T, geometry = cplanar)
