library(rgbif)
library(sf)
library(sp)
library(terra)
library(eSDM)

###############
#bring in data and get it usable

COHA = read.csv("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/cohaprey.csv")
COHAcoords = subset(COHA, select = latitude:longitude)
COHAprey = subset(COHA, select = prey_desc:prey_class)
COHA = cbind(COHAcoords, COHAprey)

#remove s, m, l avian
COHA <- subset(COHA, prey_desc != "small avian")
COHA <- subset(COHA, prey_desc != "medium avian")
COHA <- subset(COHA, prey_desc != "large avian")

COHAaves = subset(COHA, prey_class = "ave")
COHAmammals = subset(COHA, prey_class = "mammalia")
COHAherp = subset(COHA, prey_class = "sauria")

#US map
#https://hub.arcgis.com/datasets/1b02c87f62d24508970dc1a6df80c98e/explore?location=38.753686%2C-99.481986%2C3.89
us <- terra::vect("C:/Users/wgibs/OneDrive/Desktop/COHA_SSHA Diet/States_shapefile.shp")
us <- us[us$State_Code != 'AK' & us$State_Code != 'HI',] # remove AK and HI
plot(us)
crs(us)

# reproject in lambert conic projection; makes creating grids (in km) easier
lam <- "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"
us1 <- terra::project(us, crs(lam))
plot(us1, axes = T)

################
#creating a polygon around prey points

birdpoly = pts2poly_centroids(COHAaves, 25, crs = crs(us1))
birdpoly = st_as_sf(birdpoly$geometry)
birdpoly = as_Spatial(birdpoly)
mammalpoly = pts2poly_centroids(COHAmammals, 25, crs = crs(us1))
herppoly = pts2poly_centroids(COHAherp, 25, crs = crs(us1))

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
           hasCoordinate = T, geometry = samp)