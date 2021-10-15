#########################################################################################################################
# R Code Case Study 2
# Enhancing Animal Movement Analyses: Spatiotemporal Matching of Animal Positions with Remotely Sensed Data 
# Using Google Earth Engine and R   
# Ramiro D. Crego, Majaliwa M.  Masolele, Grant Connette, Jared A. Stabach
#########################################################################################################################

#Libraries
require(sp)
require(raster)
require(rgdal)
require(rgeos)
require(maptools)
require(mapview)
require(sf)
require(dplyr)
require(plyr)
library(move)
library(adehabitatLT)
library(adehabitatHR)
library(ggplot2)
# Initialize rgee
library(rgee)
ee_Initialize()


## Load Wildebeest data
wild <- read.csv("./Data/White-bearded wildebeest in Kenya.csv", header=TRUE)

# Subsetting only timestamp,lat, long and unique id column
wild <- wild[,c(3:6,13)]

#Checking the structure of the data
str(wild)

# Adding date column
wild$Date <- as.Date(wild$timestamp)

# Relocating the id to be the first column
wild <- wild %>% 
  mutate(id = tag.local.identifier) %>% 
  relocate(id)

#Summarizing the data by each id to identify the starting date and ending date
ddply(wild, ~id, summarise, min = min(Date), max = max(Date))

# Subsetting the data only from 2011-01-01 to 2012-01-01
wild <- subset(wild, Date >= "2011-01-01" & Date < "2012-01-01")

# Format the timestamp. This dataset is in UTC, and we need UTC to match with images from the ERA5_LAND product.
wild$timestamp <- as.POSIXct(wild$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC")

#Checking now if the timezone is in UTC
attr(wild$timestamp, "tzone")

# Create trajectory
temp <- as.matrix(cbind(wild$location.long, wild$location.lat))

# Project to UTM 36S
xy <- project(temp, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m")

# Create a reference date and use setNA to re-run trajectory
refda <- strptime("00:00", "%H:%M", tz="UTC")

traj.raw <-as.ltraj(xy = xy, date = wild$timestamp, id = wild$id, typeII = TRUE, slsp = c("remove"))

# Create NA values based on refda
traj.NA <- setNA(traj.raw, refda, 3, units = "hour") 

#Make a regular trajectory based on refda
traj.reg <- sett0(traj.NA, refda, 3, units = "hour")

#Checking whether the trajectory is regular or irregular
is.regular(traj.reg)

#Converting the trajectory to a dataframe
wild2<-ld(traj.reg)
ddply(wild2, ~id, summarise, min = min(date), max = max(date))

# Select animals with a full year of data

# Randomly selected 10 animals
wildebeest<- wild2 %>% filter(id %in% c('2829','2832','2836','2844','30069','30072','30075','30077','30082','30086','30076','30085')) %>% na.omit() %>% dplyr::select(c('x','y','date','id','dist'))
unique(wildebeest$id)
head(wildebeest)
#wildebeest<-wildebeest[1:1000,]
wildsf <- st_as_sf(wildebeest, coords = c('x','y'), crs="+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m")
wildsf <- st_transform(wildsf, crs = 4326)
#Visualize data
mapview(wildsf)
wildsf

#Set the date as a string with a ‘YYYY-MM-DDTHH:MM:SS’.
wildsf$date <- as.factor(wildsf$date)
wildsf$Date <- sub(" ", "T", wildsf$date)
str(wildsf)



# Load hourly climatic data
start<-"2010-12-30"
end<-"2012-01-01"
imagecoll<-ee$ImageCollection('ECMWF/ERA5_LAND/HOURLY')$filterDate(start,end)
band <- "temperature_2m" #Name of the band to use. You can change to EVI for instance when using MOD13Q1.

#Setting the GEE functions
#Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#Join Image and Points based on a maxDifference Filter within a temporal window

#Set temporal window in days for filter. This will depend on the remote sensing data used.
tempwin <- 0.1

#Set the filter
maxDiffFilter<-ee$Filter$maxDifference(
  difference=tempwin*24*60*60*1000, #days * hr * min * sec * milliseconds
  leftField= "date_millis", #Timestamp of the telemetry data
  rightField="system:time_start" #Image date
)

# Define the join.
saveBestJoin<-ee$Join$saveBest(
  matchKey="bestImage",
  measureKey="timeDiff"
)

#Function to add property with raster pixel value from the matched image
add_value<-function(feature){
  #Get the image selected by the join
  img1<-ee$Image(feature$get("bestImage"))$select(band)
  #Extract geometry from the feature
  point<-feature$geometry()
  #Get pixel value for each point at the desired spatial resolution (argument scale)
  pixel_value<-img1$sample(region=point, scale=250, tileScale = 16, dropNulls = F) 
  #Return the data containing pixel value and image date.
  feature$setMulti(list(PixelVal = pixel_value$first()$get(band), DateTimeImage = img1$get('system:index')))
}

# Function to remove image property from features
removeProperty<- function(feature) {
  #Get the properties of the data
  properties = feature$propertyNames()
  #Select all items except images
  selectProperties = properties$filter(ee$Filter$neq("item", "bestImage"))
  #Return selected features
  feature$select(selectProperties)
}



# Extract raster value
wildsf$uniq <- rep(1:100, each=1000)[1:nrow(wildsf)] 

start_time <- Sys.time()
dataoutput <- data.frame()
for(x in unique(wildsf$uniq)){
  data1 <- wildsf %>% filter(uniq == x)
  # Send sf to GEE
  data <- sf_as_ee(data1)
  # Transform day into milliseconds
  data<-data$map(add_date)
  # Apply the join.
  Data_match<-saveBestJoin$apply(data, imagecoll, maxDiffFilter)
  #Add NDVI to the data
  DataFinal<-Data_match$map(add_value)
  #Remove image property from the data
  DataFinal<-DataFinal$map(removeProperty)
  # Transform GEE object in sf
  temp<- ee_as_sf(DataFinal)
  # append
  dataoutput <- rbind(dataoutput, temp)
}
end_time <- Sys.time()
end_time - start_time
head(dataoutput)

tempwild <- dataoutput
tempwild <- st_drop_geometry(tempwild)
head(tempwild)

names(tempwild)[3] <- band
tempwild$temperature_2m <- tempwild$temperature_2m - 273.15 # convert kelvin to celcius
head(tempwild)

#saving as csv
#write.csv(tempwild,"./wildtempdata.csv")

tempwild <- read.csv("./wildtempdata.csv", header = T)
head(tempwild)


#####################
## Figure 3

tempwild$id <- as.factor(tempwild$id)
levels(tempwild$id)[levels(tempwild$id) ==  "2829"] <- "Mara-1"
levels(tempwild$id)[levels(tempwild$id) ==  "2832"] <- "Mara-2"
levels(tempwild$id)[levels(tempwild$id) ==  "2836"] <- "Mara-3"
levels(tempwild$id)[levels(tempwild$id) ==  "2844"] <- "Mara-4"
levels(tempwild$id)[levels(tempwild$id) ==  "30069"] <- "Amboseli-1"
levels(tempwild$id)[levels(tempwild$id) ==  "30072"] <- "Athi-Kaputiei-1"
levels(tempwild$id)[levels(tempwild$id) ==  "30075"] <- "Amboseli-2"  
levels(tempwild$id)[levels(tempwild$id) ==  "30077"] <- "Athi-Kaputiei-2"
levels(tempwild$id)[levels(tempwild$id) ==  "30082"] <- "Athi-Kaputiei-3"
levels(tempwild$id)[levels(tempwild$id) ==  "30086"] <- "Athi-Kaputiei-4"
levels(tempwild$id)[levels(tempwild$id) ==  "30076"] <- "Amboseli-3"
levels(tempwild$id)[levels(tempwild$id) ==  "30085"] <- "Amboseli-4"

#Change order of levels
tempwild$id <- factor(tempwild$id, levels = c('Mara-1', 'Mara-2', 'Mara-3', 'Mara-4', 'Athi-Kaputiei-1', 'Athi-Kaputiei-2', 'Athi-Kaputiei-3', 'Athi-Kaputiei-4', 'Amboseli-1', 'Amboseli-2', 'Amboseli-3', 'Amboseli-4'))

g <- ggplot(tempwild, aes(x = Temp, y = dist)) + facet_wrap(~id, scales = "free_y") +   
  coord_cartesian(ylim = c(0,2500)) +
  geom_point(alpha = 0.4, colour = "grey30") +
  ylab("Step distance (m)") + xlab("Temperature (°C)") +
  stat_smooth(aes(x = Temp, y = dist, colour = id, fill = id), method = "lm", formula = y ~ x, se = T) + 
  scale_color_manual(values = c(
    'Mara-1' = 'green4', 'Mara-2' = 'green4', 'Mara-3' = 'green4', 'Mara-4' = 'green4', 'Athi-Kaputiei-1' = 'blue', 'Athi-Kaputiei-2' = 'blue', 'Athi-Kaputiei-3' = 'blue', 'Athi-Kaputiei-4' = 'blue', 'Amboseli-1' = 'red', 'Amboseli-2' = 'red', 'Amboseli-3' = 'red', 'Amboseli-4' = 'red')) +
  scale_fill_manual(values = c(
    'Mara-1' = 'green4', 'Mara-2' = 'green4', 'Mara-3' = 'green4', 'Mara-4' = 'green4', 'Athi-Kaputiei-1' = 'blue', 'Athi-Kaputiei-2' = 'blue', 'Athi-Kaputiei-3' = 'blue', 'Athi-Kaputiei-4' = 'blue', 'Amboseli-1' = 'red', 'Amboseli-2' = 'red', 'Amboseli-3' = 'red', 'Amboseli-4' = 'red')) +
  theme_classic() +
  theme(text = element_text(size=13), legend.position = c(2, 2), strip.background.x = element_rect(color = "white"))
g


jpeg(filename = "Fig3.jpeg", width = 10, height = 6, units = "in", res  = 600, bg = "white")
g
dev.off()



### LMEM Bayesian analysis using JAGS
library(jagsUI)
library(MCMCvis)

tempwild$ind <- as.numeric(tempwild$id)

# Bundle and summarize the data set passed to JAGS

mn <- mean(tempwild$Temp)
sd <- sd(tempwild$Temp)
Temperature <- (tempwild$Temp - mn) / sd
hist(Temperature, col = "grey")

str(bdata <- list(dist = as.numeric(log(tempwild$dist)),  pop = as.numeric(tempwild$ind), temp = Temperature, ngroups = max(as.numeric(tempwild$ind)), n = nrow(tempwild)))

# Specify model in BUGS language
cat(file = "lme.model.txt", "
model {

# Priors
 for (i in 1:ngroups){		
    alpha[i] ~ dnorm(mu.int, tau.int)	# Random intercepts
    beta[i] ~ dnorm(mu.slope, tau.slope)# Random slopes
 }

 mu.int ~ dnorm(0, 0.001)		# Mean hyperparameter for random intercepts
 tau.int <- 1 / (sigma.int * sigma.int)
 sigma.int ~ dunif(0, 100)		# SD hyperparameter for random intercepts

 mu.slope ~ dnorm(0, 0.001)		# Mean hyperparameter for random slopes
 tau.slope <- 1 / (sigma.slope * sigma.slope)
 sigma.slope ~ dunif(0, 100)		# SD hyperparameter for slopes

 tau <- 1 / ( sigma * sigma)		# Residual precision
 sigma ~ dunif(0, 100)			# Residual standard deviation

# Likelihood
 for (i in 1:n) {
    dist[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha[pop[i]] + beta[pop[i]]* temp[i]
 }
}
")

# Inits function
inits <- function(){ list(alpha = rnorm(12, 0, 2), beta = rnorm(12, 0, 2), 
                          mu.int = rnorm(1, 0, 1), sigma.int = rlnorm(1), mu.slope = rnorm(1, 0, 1), 
                          sigma.slope = rlnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
params <- c("alpha", "beta", "mu.int", "sigma.int", "mu.slope", "sigma.slope", "sigma")

# MCMC settings
nc <- 3  ;  ni <- 5000  ;  nb <- 1500  ;  nt <- 5

# Call JAGS, check convergence and summarize posteriors
out <- jags(bdata, inits, params, "lme.model.txt", n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel = TRUE)
plot(out)
print(out, dig = 3)
out
