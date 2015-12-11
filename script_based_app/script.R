


# Third example: nonwhites in Glasgow only --------------------------------

# N : (uppercase) : the total population count in the area



rm(list = ls())
# Parameters --------------------------------------------------------------


# Please change this to the location of the base directory. The base directory will 
# contain a file name 'script.R', and the subdirectories 
# input_data
# output_data
# scripts

base_dir_location <- "E:/repos/dissimilarity_shiny_app/script_based_app"

# Please adjust the following model parameters if required 

n_burnin <- 2000
n_sample <- 10000

n_posterior <- 1000 # the number of estimates of D to extract


# The name and location of the file to save samples to
outfile <- "E:/repos/dissimilarity_shiny_app/script_based_app/output_data/glasgow_nonwhite_total.csv"





# Load scripts and packages  ----------------------------------------------

# Server - load scripts  --------------------------------------------------

setwd(base_dir_location)
source("scripts/LoadPackages.R")

RequiredPackages(
  required = c(
    "spdep", "rgdal", 
    "CARBayes","Rcpp", "MASS", 
    "readr", "dplyr"
  )
)              


source("scripts/D_compute.r")
sourceCpp("scripts/cppfunctions.cpp")
# 



# Load input data ---------------------------------------------------------

# Load shapefile with attributes

shp_with_attributes <- readOGR(
  dsn = "input_data/scotland_nonwhite_total",
  layer = "scotland_dz_nonwhite_total"
)

# Load dz to local authority lookup 

dzla_lookup <- read.csv("input_data/la_to_dz.csv") 

dzs_in_glasgow_city <- subset(dzla_lookup, local_authority == "Glasgow City")$datazone 
dzs_in_glasgow_city <- as.character(dzs_in_glasgow_city)
  
shp_ss <- subset(shp_with_attributes, subset = geo_id %in% dzs_in_glasgow_city)

# Now to calculate neighbourhood matrix from the above

w_nb <- poly2nb(shp_ss, row.names = shp_ss@data$geo_id)


tmp <- card(w_nb) #shows the number of neighbours for each areal unit
names(tmp) = shp_ss@data$geo_id
geo_id_without_neighbours <- names(tmp)[which(tmp == 0)]
# None! 


w_nb <- nb2mat(w_nb, style = "B", zero.policy = TRUE)
# 2) subset w_nb to exclude areas with no neigbours
dim(w_nb)

w_nb <- w_nb[which(tmp > 0), which(tmp > 0)]
dim(w_nb)



mdl <- S.CARiar(
  formula= numer  ~ 1,
  trials = shp_ss$total,
  W=w_nb,
  data=shp_ss,
  family="binomial",
  burnin = n_burnin,
  n.sample = n_sample, 
  verbose = TRUE
)

out <- list(
  datazones = dta_dz,
  denominator = dta$denominator,
  beta = mdl$samples$beta[,1],
  phi = mdl$samples$phi
)

# This time, it works!

K <- n_posterior      
phi <- mdl$samples$phi
beta <- mdl$samples$beta[,1]
denominator <- shp_ss$total

out <- array(NA, K)
for(k in 1:K){
  p.current <- exp(phi[k ,] + beta[k])   / (1 + exp(phi[k ,] + beta[k]))
  p.current.overall <- sum(p.current * denominator) / sum(denominator)
  out[k] <- sum(denominator * abs(p.current - p.current.overall)) / 
    ( 2 * sum(denominator) * p.current.overall * (1-p.current.overall))                 
}


# Out contains (initially) 1000 estimates of the true value of D

write.csv(out, file = outfile, row.names = F)


# Examples for the whole of Scotland, which  fail at the model stage -----------------------------------


# 
# # 10/12/2015
# 
# # App to produce posterior values for dissimilarity scores based on Lee's model
# 
# 
# # Instructions to users: 
# 
# # please put shapefiles in separate folders within 
# # input_data/shapefiles
# 
# # Please ensure these shapefile folders include the necessary files (.shp, .shx, .dbf)
# # IMPORTANT:
# # Please ensure the .dbf file includes two variables:
# # numer :  the numerator count 
# # total :  the total population count in the area
# 
# 
# 
# rm(list = ls())
# # Parameters --------------------------------------------------------------
# 
# 
# # Please change this to the location of the base directory. The base directory will 
# # contain a file name 'script.R', and the subdirectories 
# # input_data
# # output_data
# # scripts
# 
# base_dir_location <- "E:/repos/dissimilarity_shiny_app/script_based_app"
# 
# # Please adjust the following model parameters if required 
# 
# n_burnin <- 2000
# n_sample <- 10000
# n_posterior <- 1000
# 
# 
# 
# # Load scripts and packages  ----------------------------------------------
# 
# # Server - load scripts  --------------------------------------------------
# 
# setwd(base_dir_location)
# source("scripts/LoadPackages.R")
# 
# RequiredPackages(
#   required = c(
#     "spdep", "rgdal", 
#     "CARBayes","Rcpp", "MASS", 
#     "readr", "dplyr"
#   )
# )              
# 
# 
# source("scripts/D_compute.r")
# sourceCpp("scripts/cppfunctions.cpp")
# # 
# 
# 
# 
# # Load input data ---------------------------------------------------------
# 
# # Load shapefile with attributes
# 
# shp_with_attributes <- readOGR(
#   dsn = "input_data/scotland_noreligion_total",
#   layer = "scotland_dz_n_noreligion_N_total"
# )
# 
# 
# 
# # Now to calculate neighbourhood matrix from the above
# 
# w_nb <- poly2nb(shp_with_attributes, row.names = shp_with_attributes@data$geo_id)
# 
# 
# tmp <- card(w_nb) #shows the number of neighbours for each areal unit
# names(tmp) = shp_with_attributes@data$geo_id
# geo_id_without_neighbours <- names(tmp)[which(tmp == 0)]
# 
# # to do: 
# # 1) drop geo_ids in shp_with_attributes where there are no neighbours 
# shp_ss <- subset(shp_with_attributes, !(geo_id %in% geo_id_without_neighbours)) 
# 
# # This drops the length to 6492
# 
# # 2) subset w_nb to exclude areas with no neigbours
# 
# w_nb <- nb2mat(w_nb, style = "B", zero.policy = TRUE)
# # 2) subset w_nb to exclude areas with no neigbours
# dim(w_nb)
# 
# w_nb <- w_nb[which(tmp > 0), which(tmp > 0)]
# dim(w_nb)
# 
# 
# 
# mdl <- S.CARiar(
#   formula= numer  ~ 1,
#   trials = shp_ss$total,
#   W=w_nb,
#   data=shp_ss,
#   family="binomial",
#   burnin = n_burnin,
#   n.sample = n_sample, 
#   verbose = TRUE
# )
# 
# # # This produces the following error: 
# # Setting up the model
# # Collecting 10000 samples
# # |==========                                                                               |  11%
# # Error in if (prob > runif(1)) { : missing value where TRUE/FALSE needed
# 
# # Which seems to be quite a low level error - buried within a function within a function
# # within CARBayes
# 
# # I will try another example 
# 
# 
# # Second example: nonwhite/total  -----------------------------------------
# 
# 
# rm(list = ls())
# # Parameters --------------------------------------------------------------
# 
# 
# # Please change this to the location of the base directory. The base directory will 
# # contain a file name 'script.R', and the subdirectories 
# # input_data
# # output_data
# # scripts
# 
# base_dir_location <- "E:/repos/dissimilarity_shiny_app/script_based_app"
# 
# # Please adjust the following model parameters if required 
# 
# n_burnin <- 2000
# n_sample <- 10000
# n_posterior <- 1000
# 
# 
# 
# # Load scripts and packages  ----------------------------------------------
# 
# # Server - load scripts  --------------------------------------------------
# 
# setwd(base_dir_location)
# source("scripts/LoadPackages.R")
# 
# RequiredPackages(
#   required = c(
#     "spdep", "rgdal", 
#     "CARBayes","Rcpp", "MASS", 
#     "readr", "dplyr"
#   )
# )              
# 
# 
# source("scripts/D_compute.r")
# sourceCpp("scripts/cppfunctions.cpp")
# # 
# 
# 
# 
# # Load input data ---------------------------------------------------------
# 
# # Load shapefile with attributes
# 
# shp_with_attributes <- readOGR(
#   dsn = "input_data/scotland_nonwhite_total",
#   layer = "scotland_dz_nonwhite_total"
# )
# 
# 
# 
# # Now to calculate neighbourhood matrix from the above
# 
# w_nb <- poly2nb(shp_with_attributes, row.names = shp_with_attributes@data$geo_id)
# 
# 
# tmp <- card(w_nb) #shows the number of neighbours for each areal unit
# names(tmp) = shp_with_attributes@data$geo_id
# geo_id_without_neighbours <- names(tmp)[which(tmp == 0)]
# 
# # to do: 
# # 1) drop geo_ids in shp_with_attributes where there are no neighbours 
# shp_ss <- subset(shp_with_attributes, !(geo_id %in% geo_id_without_neighbours)) 
# 
# # This drops the length to 6492
# 
# # 2) subset w_nb to exclude areas with no neigbours
# 
# w_nb <- nb2mat(w_nb, style = "B", zero.policy = TRUE)
# # 2) subset w_nb to exclude areas with no neigbours
# dim(w_nb)
# 
# w_nb <- w_nb[which(tmp > 0), which(tmp > 0)]
# dim(w_nb)
# 
# 
# 
# mdl <- S.CARiar(
#   formula= numer  ~ 1,
#   trials = shp_ss$total,
#   W=w_nb,
#   data=shp_ss,
#   family="binomial",
#   burnin = n_burnin,
#   n.sample = n_sample, 
#   verbose = TRUE
# )
# 
# # This also breaks: 
# # Setting up the model
# # Collecting 10000 samples
# # |====================                                                                     |  22%
# # Error in if (prob > runif(1)) { : missing value where TRUE/FALSE needed


