##### R-Code for Poster: Population grids protected by the Cell Key Method ### #
#
# 01 - Download and prepare source data
#
# contact: Martin Möhler, martin.moehler@destatis.de
# 2024-04-09
############################################################################## #


library(data.table)
library(sf)
library(dplyr) 


# customize if needed:
data_folder <- file.path(getwd(), "R/AnoSiDat2024") 


# ----- get data from web -----

# Zensus 2011 population by 1ha grid cell 
web_pop <- paste0("https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/",
                  "DemografischeGrunddaten/csv_Bevoelkerung_100m_Gitter.zip",
                  "?__blob=publicationFile&v=2")

# CORINE Land Cover 5 ha, Stand 2018 (CLC5-2018)
web_clc <- paste0("https://daten.gdz.bkg.bund.de/produkte/dlm/clc5_2018/aktuell/",
                  "clc5_2018.utm32s.shape.zip")

# download
download.file(web_pop, paste0(data_folder, "/csv_Bevoelkerung_100m_Gitter.zip"),
              mode = "wb") # ~ 105 MB
download.file(web_clc, paste0(data_folder, "/clc5_2018.utm32s.shape.zip"),
              mode = "wb") # ~ 1.3 GB

# unzip
unzip(paste0(data_folder, "/csv_Bevoelkerung_100m_Gitter.zip"), exdir = data_folder)
unzip(paste0(data_folder, "/clc5_2018.utm32s.shape.zip"), exdir = data_folder)


# ----- read in and crop  data -----

## population grid data (csv)

pop <- fread(file.path(data_folder, "Zensus_Bevoelkerung_100m-Gitter.csv")) 

pop <- pop[pop$x_mp_100m > 4100000 & pop$x_mp_100m < 4200000 &
             pop$y_mp_100m > 3100000 & pop$y_mp_100m < 3200000, ]

pop$Einwohner[pop$Einwohner == - 1] <- 0

# expand to record-level microdata
pop_grid <- pop[pop$Einwohner != 0, 
                .(1:Einwohner, x_mp_100m, y_mp_100m), 
                by = .(Gitter_ID_100m)] %>% select(-V1)
pop_grid$person_id <- 1:nrow(pop_grid)


## land cover data (shape files)

shp <- vector("list", 5)
for(i in 1:5) {
  
  clc5_name <- paste0("clc5/clc5_class", i, "xx.shp")
  
  shp[[i]] <- st_read(file.path(data_folder, "/clc5_2018.utm32s.shape/", clc5_name)) %>%
    st_transform(crs = "EPSG:3035") %>%
    st_crop(xmin = min(pop$x_mp_100m) - 50, xmax = max(pop$x_mp_100m) + 50,
            ymin = min(pop$y_mp_100m) - 50, ymax = max(pop$y_mp_100m) + 50)
}

clc18 <- bind_rows(shp)[, "CLC18"]


# ----- store data -----

save(pop, pop_grid, file = file.path(data_folder, "/pop_data_cell29.RData"), 
     compress = TRUE)

st_write(clc18, file.path(data_folder, "clc18.shp"))

rm(list = ls())

