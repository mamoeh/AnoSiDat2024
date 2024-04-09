##### R-Code for Poster: Population grids protected by the Cell Key Method ### #
#
# 02 - Apply protection and analyse differences
#
# contact: Martin Möhler, martin.moehler@destatis.de
# 2024-04-09
############################################################################## #


library(ptable)
library(raster)
library(dplyr)
library(sdcSpatial)
library(ggplot2)


source("R/AnoSiDat2024/functions_ckm_raster.R") # helper functions

load("R/AnoSiDat2024/pop_data_cell29.RData") # grid data
clc18 <- st_read("R/AnoSiDat2024/clc18.shp") # land cover data

cbpal <- c("#0072B2", "#E69F00", "#56B4E9", "#D55E00") # color palette


# ----- 1. prepare for CKM -----

# make p-table
# see: https://cran.r-project.org/web/packages/ptable/vignettes/introduction.html
pt1 <- create_cnt_ptable(D = 4, V = 1.8, js = 1)
pt2 <- create_cnt_ptable(D = 6, V = 2.5, js = 2)

pt1@pTable$variant <- "I"
pt2@pTable$variant <- "II"


# visualize ptables together
ptabs <- pt1@pTable[pt1@pTable$i > 0, ]
for(i in max(pt1@pTable$i):max(pt2@pTable$i)) {
  pt    <- pt1@pTable[pt1@pTable$i == max(pt1@pTable$i), ]
  pt$i  <- i
  pt$j  <- pt$i + pt$v
  ptabs <- rbind(ptabs, pt)
}
ptabs <- rbind(ptabs, pt2@pTable[pt2@pTable$i > 0, ])
ptabs$i <- factor(ptabs$i, levels = 1:max(ptabs$i), 
                  labels = paste("i:", c(1:(max(ptabs$i) - 1), 
                                         paste0(max(ptabs$i), "+"))))

ggplot(ptabs, 
       aes(x = as.integer(v), y = p, group = variant)) +
  geom_point(aes(fill = variant), pch = 21, alpha = 0.5, cex = 1) +
  geom_line(aes(color = variant), alpha = 0.5) +
  scale_fill_manual(values = c(cbpal[1], cbpal[2])) +
  scale_color_manual(values = c(cbpal[3], cbpal[4])) +
  facet_wrap(~i) +
  theme_bw() +
  xlab("noise value added") +
  ylab("probability") +
  theme(strip.background = element_rect(fill = cbpal[3]),
        legend.position  = "bottom")


# ----- 2. apply CKM -----

# draw record keys
set.seed(20240307)
pop_grid$rk <- runif(nrow(pop_grid))

# aggregate microdata to raster
r_dat <- sdc_raster_ckm(pop_grid[, c("x_mp_100m", "y_mp_100m")], 
                        r = raster(resolution = 100, 
                                   xmn = min(pop_grid$x_mp_100m) - 50,
                                   xmx = max(pop_grid$x_mp_100m) + 50,
                                   ymn = min(pop_grid$y_mp_100m) - 50,
                                   ymx = max(pop_grid$y_mp_100m) + 50), 
                        variable = 1, min_count = 3,
                        rkey = pop_grid$rk)

# perturbation step
r_ckm1 <- protect_ckm(r_dat, ptab = pt1)
r_ckm2 <- protect_ckm(r_dat, ptab = pt2)

# select smaller focus regions
fa1 <- extent(r_dat$value, 425, 475, 25, 75)
fa2 <- extent(r_dat$value, 210, 260, 670, 720)
fa3 <- extent(r_dat$value, 925, 975, 475, 525)


## inspect map

r_df <- as.data.frame(r_dat$value$count, xy = TRUE) %>% na.omit()

ggplot(r_df, aes(x, y, fill = count)) +
  geom_raster() +
  scale_fill_viridis_c(direction =-1) +
  theme_minimal() +
  ggtitle("LAEA 100kmN31E41") +
  xlab("E") +
  ylab("N") +
  geom_rect(aes(xmin = fa1@xmin, 
                xmax = fa1@xmax, 
                ymin = fa1@ymin, 
                ymax = fa1@ymax), fill = NA, color = "darkred") +
  geom_rect(aes(xmin = fa2@xmin, 
                xmax = fa2@xmax, 
                ymin = fa2@ymin, 
                ymax = fa2@ymax), fill = NA, color = "darkred") +
  geom_rect(aes(xmin = fa3@xmin, 
                xmax = fa3@xmax, 
                ymin = fa3@ymin, 
                ymax = fa3@ymax), fill = NA, color = "darkred")


###### Analyses ########################################################

# ----- Analysis 1: Landcover types -----

r_ckm_fa <- rbind(as.data.frame(crop(r_dat$value$count, fa1), xy = TRUE) %>%
                    mutate(variant = "orig.", fa = "area 1"),
                  as.data.frame(crop(r_dat$value$count, fa2), xy = TRUE) %>%
                    mutate(variant = "orig.", fa = "area 2"),
                  as.data.frame(crop(r_dat$value$count, fa3), xy = TRUE) %>%
                    mutate(variant = "orig.", fa = "area 3"),
                  as.data.frame(crop(r_ckm1$value$count, fa1), xy = TRUE) %>%
                    mutate(variant = "I", fa = "area 1"),
                  as.data.frame(crop(r_ckm1$value$count, fa2), xy = TRUE) %>%
                    mutate(variant = "I", fa = "area 2"),
                  as.data.frame(crop(r_ckm1$value$count, fa3), xy = TRUE) %>%
                    mutate(variant = "I", fa = "area 3"),
                  as.data.frame(crop(r_ckm2$value$count, fa1), xy = TRUE) %>%
                    mutate(variant = "II", fa = "area 1"),
                  as.data.frame(crop(r_ckm2$value$count, fa2), xy = TRUE) %>%
                    mutate(variant = "II", fa = "area 2"),
                  as.data.frame(crop(r_ckm2$value$count, fa3), xy = TRUE) %>%
                    mutate(variant = "II", fa = "area 3")) %>%
  na.omit() %>%
  st_as_sf(coords = c("x", "y"))

st_crs(r_ckm_fa) <- "EPSG:3035"

r_ckm_fa <- st_join(r_ckm_fa, clc18, st_intersects)

r_ckm_fa$landcover <- 4
r_ckm_fa$landcover[substr(r_ckm_fa$CLC18, 1, 1) %in% c("2", "3")] <- 3
r_ckm_fa$landcover[r_ckm_fa$CLC18 == "111"] <- 1
r_ckm_fa$landcover[r_ckm_fa$CLC18 == "112"] <- 2
r_ckm_fa$landcover <- factor(r_ckm_fa$landcover, levels = 1:4, 
                             labels = c("cU", "dcU", "ff", "oth"))
r_ckm_fa$variant <- factor(r_ckm_fa$variant, levels = c("orig.", "I", "II"))

lc_aggr <- st_drop_geometry(r_ckm_fa) %>% group_by(variant, fa, landcover) %>%
  summarise(count = sum(count))

ggplot(lc_aggr, aes(x = landcover, y = count, fill = variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("grey", cbpal[1], cbpal[2])) +
  facet_wrap(~fa) +
  xlab(NULL) +
  theme_bw() +
  theme(strip.background   = element_rect(fill = cbpal[3]),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())


# ----- Analysis 2: Buffer operations -----

buff_rad <- 1000 # buffer radius in meters

# exclude raster margin for picking buffer centers
nr <- buff_rad / 100
rcs <- c(1:nr, (ncol(r_dat$value) - nr):ncol(r_dat$value))

r_dat$value$inner <- TRUE
r_dat$value$inner[c(cellFromCol(r_dat$value, rcs), cellFromRow(r_dat$value, rcs))] <- FALSE

# raster to coordinates
buff_coords       <- as.data.frame(r_dat$value$inner, xy = TRUE)
buff_coords$count <- getValues(r_dat$value$count)
buff_coords$ckm1  <- getValues(r_ckm1$value$count)
buff_coords$ckm2  <- getValues(r_ckm2$value$count)
buff_coords$prob  <- buff_coords$count / sum(buff_coords$count, na.rm = TRUE)

# pick buffer centers

nSamp <- 100 # sample size
cands <- which(!is.na(buff_coords$count)) # cells eligible for buffer center

set.seed(20240328)
bcents <- sample(cands, nSamp, replace = FALSE)

# (empty) results data frame
buff_results <- data.frame(i = rep(1:nSamp, 2),
                           x = rep(buff_coords$x[bcents], 2), 
                           y = rep(buff_coords$y[bcents], 2),
                           variant = rep(c("I", "II"), each = nSamp),
                           count_true = 0, count_ckm = 0)

for(i in 1:nSamp) {
  
  buff_coords$dist <- 
    sqrt((buff_results$x[i] - buff_coords$x)^2 + (buff_results$y[i] - buff_coords$y)^2)
  buff_coords$in_buff <- buff_coords$dist <= buff_rad
  
  buff_results$count_true[buff_results$i == i] <- 
    sum(buff_coords$count[buff_coords$in_buff], na.rm = TRUE)
  buff_results$count_ckm[buff_results$i == i & buff_results$variant == "I"] <- 
    sum(buff_coords$ckm1[buff_coords$in_buff], na.rm = TRUE)
  buff_results$count_ckm[buff_results$i == i & buff_results$variant == "II"] <- 
    sum(buff_coords$ckm2[buff_coords$in_buff], na.rm = TRUE)
}

buff_results$diff_abs <- 
  abs(buff_results$count_ckm - buff_results$count_true) / buff_results$count_true * 100

ggplot(buff_results[c(101:200, 1:100), ], aes(count_true, count_ckm, color = variant)) +
  geom_abline(intercept = c(0, 0), slope = 1, color = "grey50") +
  geom_point(cex = .5) +
  scale_color_manual(values = c(cbpal[1], cbpal[2])) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  xlab("count without CKM") +
  ylab("count with CKM") +
  theme_bw()

ggplot(buff_results, aes(count_true, diff_abs, color = variant)) +
  geom_point(cex = .5, show.legend = FALSE) +
  scale_color_manual(values = c(cbpal[1], cbpal[2])) +
  scale_x_continuous(trans = "log10") +
  xlab("count without CKM") +
  ylab("relative error (%)") +
  theme_bw()


# ----- Analysis 3: Hot spots -----

## helper functions

find_hotspots <- function(x, q) {
  
  threshold <- quantile(getValues(x), q, na.rm = TRUE)
  hs <- x > threshold
  hs
}

jaccard_sim <- function(x, y) {
  
  xv <- getValues(x)
  yv <- getValues(y)
  
  sum(xv & yv, na.rm = TRUE) / sum(xv | yv, na.rm = TRUE)
}

hotspot_value_sum <- function(x, q) {
  
  xh <- find_hotspots(x, q)
  sum(getValues(x * xh), na.rm = TRUE)
}


## application

q <- 0.8 # quantile threshold for hot spot definition

r_dat$value$hot  <- find_hotspots(r_dat$value$count,  q)
r_ckm1$value$hot <- find_hotspots(r_ckm1$value$count, q)
r_ckm2$value$hot <- find_hotspots(r_ckm2$value$count, q)

## Jaccard similarity

# full map
jaccard_sim(r_dat$value$hot, r_ckm1$value$hot)
jaccard_sim(r_dat$value$hot, r_ckm2$value$hot)
# focus area 1
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa1), q),
            find_hotspots(crop(r_ckm1$value$count, fa1), q))
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa1), q),
            find_hotspots(crop(r_ckm2$value$count, fa1), q))
# focus area 2
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa2), q),
            find_hotspots(crop(r_ckm1$value$count, fa2), q))
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa2), q),
            find_hotspots(crop(r_ckm2$value$count, fa2), q))
# focus area 3
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa3), q),
            find_hotspots(crop(r_ckm1$value$count, fa3), q))
jaccard_sim(find_hotspots(crop(r_dat$value$count,  fa3), q),
            find_hotspots(crop(r_ckm2$value$count, fa3), q))


## population in hotspot

# full map
pop_hs_dat  <- hotspot_value_sum(r_dat$value$count,  q)
pop_hs_ckm1 <- hotspot_value_sum(r_ckm1$value$count, q)
pop_hs_ckm2 <- hotspot_value_sum(r_ckm2$value$count, q)
(pop_hs_ckm1 - pop_hs_dat) / pop_hs_dat * 100
(pop_hs_ckm2 - pop_hs_dat) / pop_hs_dat * 100
# focus area 1
pop_hs_dat  <- hotspot_value_sum(crop(r_dat$value$count,  fa1),  q)
pop_hs_ckm1 <- hotspot_value_sum(crop(r_ckm1$value$count, fa1),  q)
pop_hs_ckm2 <- hotspot_value_sum(crop(r_ckm2$value$count, fa1),  q)
(pop_hs_ckm1 - pop_hs_dat) / pop_hs_dat * 100
(pop_hs_ckm2 - pop_hs_dat) / pop_hs_dat * 100
# focus area 2
pop_hs_dat  <- hotspot_value_sum(crop(r_dat$value$count,  fa2),  q)
pop_hs_ckm1 <- hotspot_value_sum(crop(r_ckm1$value$count, fa2),  q)
pop_hs_ckm2 <- hotspot_value_sum(crop(r_ckm2$value$count, fa2),  q)
(pop_hs_ckm1 - pop_hs_dat) / pop_hs_dat * 100
(pop_hs_ckm2 - pop_hs_dat) / pop_hs_dat * 100
# focus area 3
pop_hs_dat  <- hotspot_value_sum(crop(r_dat$value$count,  fa3),  q)
pop_hs_ckm1 <- hotspot_value_sum(crop(r_ckm1$value$count, fa3),  q)
pop_hs_ckm2 <- hotspot_value_sum(crop(r_ckm2$value$count, fa3),  q)
(pop_hs_ckm1 - pop_hs_dat) / pop_hs_dat * 100
(pop_hs_ckm2 - pop_hs_dat) / pop_hs_dat * 100

