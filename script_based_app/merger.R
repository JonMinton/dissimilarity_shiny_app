rm(list = ls())
setwd("E:/repos/dissimilarity_shiny_app/")

require(shapefiles)
require(rgdal)
require(readr)
require(tidyr)
require(dplyr)
require(CARBayes)

shapefile <- read.shp(
  shp.name = "single_file_app/shp/scotland_2001_datazones/scotland_dz_2001.shp"
)

dtabase <- read.dbf(
  dbf.name = "single_file_app/shp/scotland_2001_datazones/scotland_dz_2001.dbf"
)

dtabase$dbf <- dtabase$dbf[,c(2,1,3,4,5)]
# moves zonecode (i.e. datazone, to first column in dtabase$dbf)


attributes <- read_csv("single_file_app/data/rel_2001.csv")

attributes <- attributes %>% 
  dplyr::select(-year)   %>% 
  spread(key = type, value = count)  %>% 
  dplyr::select(datazone, denominator = all_people, numerator = none)


class(attributes) <- "data.frame"
attributes <- attributes[!is.na(attributes$datazone),]

row.names(attributes) <- attributes$datazone

joined_data <- combine.data.shapefile(
  data = attributes, shp = shapefile, dbf = dtabase
)



writeOGR(joined_data, dsn = "script_based_app/input_data/scotland_noreligion_total", 
         layer = "scotland_dz_n_noreligion_N_total",
         driver = "ESRI Shapefile",
         morphToESRI = TRUE
         )


