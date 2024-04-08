library(sp); library(raster); library(rgeos); library(geosphere)

# read in files
coha_all <- read.csv("./COHA_RawData.csv") # read in COHA raw data from working directory
prey_table <- read.csv("./PreyWeightTable_21jan2024.csv") # read in COHA prey weights from working directory
prey_table <- prey_table[ ,c(1,3,4,6)] # subset to get columns 1, 3, 4, and 6 (prey ID, min/max weight, and avian size class)

# histogram of prey item weights
prey_tableDATA <- subset(prey_table, min_weight != "na" & min_weight != 0)
prey_tableDATA$min_weight <- as.numeric(prey_tableDATA$min_weight)
prey_tableDATA$max_weight <- as.numeric(prey_tableDATA$max_weight)
prey_tableDATA$weight <- 0  # blank column to hold randomly assigned prey weights

for(j in 1:nrow(prey_tableDATA)){
  # j = 1
  prey_j <- prey_tableDATA[j,]
  low1 <- prey_j$min_weight # lower bound of prey distribution
  high1 <- prey_j$max_weight # upper bound of prey distribution
  mean1 <- mean(c(low1, high1)) # mean is the midpoint
  samp <- round(rnorm(n = 1, mean = mean1, sd = (high1-mean1)/2),0) # sample one point, round
  samp <- ifelse(samp < 0, mean1, samp)
  prey_j$weight <- samp # add new random weight to prey
  prey_tableDATA[j,] <- prey_j # add column with new weight back to original data.frame
  print(paste0("row ", j, " done. Prey ", prey_j$prey_desc, " assigned weight of ", samp, " g"))
}

# plotting
plot(-1,-1, ylim = c(0, 0.006), xlim = c(0,500), xlab = "Prey Weight", ylab = "Probability")
c1 <- subset(prey_tableDATA, weight < 500) # object "c1" is just for the histogram on the next line
lines(density(c1$weight), lwd = 2, col = "black") # density plot line
rug(c1$weigh)

