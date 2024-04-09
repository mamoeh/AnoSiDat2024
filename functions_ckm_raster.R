# ----- Functions to integrate the Cell Key Method (CKM) with sdcSpatial -----

# Martin Möhler (martin.moehler@destatis.de)
# 2024-04-09


###### Look up the perturbed value by cell key from a count pTable ########### # 
# x_ck : length-2 vector c(cell count, cell key) for one table cell
# pt   : pTable from a ptable object
############################################################################## #

get_perturbation_value <- function(x_ck, pt) {
  
  if(is.na(x_ck[1])) { 
    NA 
  } else {
    max_x <- max(pt$i)
    pt$v[pt$i == min(max_x, x_ck[1]) & data.table::between(x_ck[2], pt$p_int_lb, pt$p_int_ub)]
  }
}


##### Create an sdc_raster object that includes cell keys #################### #
# x        : SpatialPointsDataFrame, sf (point data), or 2-column matrix of 
#            coordinates
# variable : name of field in sp or sf; when using coordinates a numeric of 
#            length x
# rkey     : name of field in sp or sf that contains record keys; when using 
#            coordinates a numeric vector of record keys with length x
############################################################################## #

sdc_raster_ckm <- function(x, variable, rkey, ...) {
  
  # default sdc_raster call
  sr <- sdcSpatial::sdc_raster(x, variable, ...)
  
  # cell keys from record keys
  if(is.numeric(rkey)) { 
    le_rkey <- nchar(rkey[1]) - 2 
  } else { le_rkey <- nchar(x[[rkey]][1]) - 2 }
  
  ck <- raster::rasterize(x, sr$value, fun = sum, field = rkey)
  ck <- round(ck %% 1, le_rkey)
  sr$value$ckeys <- ck
  sr
}


##### Perturb cell values in a raster map using CKM ########################## #
# x    : an sdc_raster object
# ptab : a ptable object
# -- currently set up to work only on count values --
############################################################################## #

protect_ckm <- function(x, ptab) {
  
  r <- x$value
  
  # extract values
  x_ck <- cbind(raster::getValues(r[["count"]]),
                raster::getValues(r[["ckeys"]]))
  
  # look up noise by cell key
  v <- apply(x_ck, 1, get_perturbation_value, pt = ptab@pTable)
  # perturb values with noise
  raster::values(r[["count"]]) <- raster::values(r[["count"]]) + v
  
  x$value <- r
  x
}

