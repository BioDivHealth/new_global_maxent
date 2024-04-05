# Load libraries -----

require(raster)
#require(rgdal)  # deprecated
require(dismo)
require(data.table)
require(sp)
#require(maptools) # deprecated
require(doParallel)
require(fasterize)
require(sf)

data(wrld_simpl) # might not work if maptools deprecated

##set raster template
template <-
 raster(
    nrow = 3600,
    ncol = 8640,
    ext = extent(-180, 180,-60, 90),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  )
#values(template) <- 1:ncell(template)

#ws_ras <- fasterize(st_as_sf(wrld_simpl), template, field = "UN")

# load pre-existing raster
ws_ras <- raster("data/wrld_simpl/wrld_simpl_raster.tif")

##mask template by area
template <- mask(template, ws_ras)

source("scripts/functions/functions6.r")

##data.table function
cbind2 <-
  function(...)
    (setattr(do.call(c, list(...)), "class", c("data.table", "data.frame")))

years = c(2030, 2050, 2070, 2080, 2010)
RCPS <- c(2.6, 8.5, 6.0, 4.5)
lcs <-
  c("gcrop",
    "gsecd",
    "gpast",
    "gurbn",
    "gothr",
    "gfsh1",
    "gfsh2",
    "gfsh3")

wrld_simpl <- read_sf("data/wrld_simpl/wrld_simpl3.shp","wrld_simpl3")
data(wrld_simpl)
wrld_simpl[wrld_simpl$NAME == "Russia", "SUBREGION"] <- 161

###read disease data
d1 <- read.csv("data/disease_table32c.csv",stringsAsFactors=FALSE)
d1$name2 <- paste(" ", d1$name, sep = "")
d1$name2 <-
  gsub(" angiostrongylus costaricensis ",
       " angiostrongylus costaricensis",
       d1$name2)
d1$spillover_rate2 <-
  cut(
    d1$cases_per_year,
    breaks = c(0, 0.99, 99, 99999, 999999999999),
    labels = FALSE
  ) / 4

###look up all diseases
load(file = "scripts/functions/dis1.r")

##read in diseases
diseases <-
  read.csv(file = "data/diseases1.csv",
           stringsAsFactors = FALSE,
           encoding = "latin1")

##get all points data for each
points1a <-
  list.files(
    "data/points",
    pattern = "all_points",
    full.names = TRUE
  )

points1 <-
  list.files(
    "data/points",
    pattern = "_points",
    full.names = TRUE
  )

points2b <- points1[!points1 %in% points1a]

points1a <-
  list.files(
    "data/points",
    pattern = "all_points",
    full.names = FALSE
  )

points1 <-
  list.files(
    "data/points",
    pattern = "_points",
    full.names = FALSE
  )

points2 <- points1[!points1 %in% points1a]
points2 <- gsub("_points.r", "", points2, fixed = TRUE)
points3 <- data.frame(species = points2, filen = points2b)
rm(points1a, points1, points2b, points2)

#load gridded livestock of the world
# TO DO: ADD YOUR OWN FILEPATH
lt2 <-
  fread("≈livestock_future_2030_2050_2070_2080b.csv")
setkey(lt2, cell_by_year)


####Get all future climate niches
# TO DO: ADD YOUR OWN FILEPATH

link2b <-
  list.files(
    "/Volumes/OS/Users/Public/Documents/resultsY2",
    pattern = ".tif",
    full.names = TRUE
  )
link2 <-
  list.files(
    "/Volumes/OS/Users/Public/Documents/resultsY2",
    pattern = ".tif",
    full.names = FALSE
  )

link2 <- gsub("_6.tif", "_6.0.tif", link2, fixed = TRUE)
link2 <-
  gsub("_present_", ";2010;present", link2, fixed = TRUE) ##check when in
link2 <- gsub("_XXX.tif", "", link2, fixed = TRUE)
link2 <- gsub(".tif", "", link2, fixed = TRUE)
link2 <- gsub("_XXX_", ";", link2, fixed = TRUE)
link2 <- gsub("_2.6", ";2.6", link2, fixed = TRUE)
link2 <- gsub("_4.5", ";4.5", link2, fixed = TRUE)
link2 <- gsub("_6.0", ";6.0", link2, fixed = TRUE)
link2 <- gsub("_8.5", ";8.5", link2, fixed = TRUE)
link2 <- gsub("_2030", ";2030", link2, fixed = TRUE)
link2 <- gsub("_2050", ";2050", link2, fixed = TRUE)
link2 <- gsub("_2070", ";2070", link2, fixed = TRUE)
link2 <- gsub("_2080", ";2080", link2, fixed = TRUE)
link2 <- gsub("-astOJIRBNlQ0", "", link2, fixed = TRUE)
link2 <- gsub("_points.r", "", link2, fixed = TRUE)
link2 <- gsub("presentXXX", "present", link2, fixed = TRUE)

link3 <- read.table(text = link2,
                    sep = ";",
                    stringsAsFactors = FALSE)

link3$filen <- link2b
names(link3)[1:3] <- c("species", "year", "RCP")
link3$RCP <- as.numeric(gsub("present", "999", link3$RCP))

###make other species regions available
link3$speciesX <- link3$species
link3$speciesX <- gsub("_1_subs_", ";", link3$speciesX, fixed = TRUE)
link3$speciesX <- gsub("_2_subs_", ";", link3$speciesX, fixed = TRUE)
link3$speciesX <- gsub("_3_subs_", ";", link3$speciesX, fixed = TRUE)
link3$speciesX <- gsub("_4_subs_", ";", link3$speciesX, fixed = TRUE)

others1 <-
  read.table(text = link3$speciesX,
             sep = ";",
             stringsAsFactors = FALSE)
names(others1) <- c("species2", "region")
link3 <- cbind(link3, others1)
rm(link2, link2b)

###
futlg <-
  read.csv(file = "data/futlg.csv", stringsAsFactors =
             FALSE)
head(futlg)
# TO DO: CHANGE THIS TO YOUR OWN FILEPATH
futlg$filen <-
  gsub("C:/", "OS/Users/Public/Documents/", futlg$filen, fixed = TRUE)
futlg$filen <- gsub("/updated_states", "", futlg$filen, fixed = TRUE)
futlg$filen <- gsub("LUH1//", "LUH1/", futlg$filen, fixed = TRUE)

## forest not forest land
fnf <- raster("data/fnf_map.txt")
plot(fnf)

##ecological data about species
spec_hab <-
  read.csv("C:\\Users\\xxxx\\Dropbox\\legion2\\habitat_both_fin.csv",
           stringsAsFactors = FALSE)
spec_hab$Snow.and.ice[spec_hab$Snow.and.ice == 0] <- NA
spec_hab$f <-
  (
    spec_hab$Evergreen.Needleleaf.forest + spec_hab$Evergreen.Broadleaf.forest +
      spec_hab$Deciduous.Needleleaf.forest + spec_hab$Deciduous.Broadleaf.forest +
      spec_hab$Mixed.forest
  ) / 5
spec_hab$nf <-
  (
    spec_hab$Closed.shrublands + spec_hab$Open.shrublands + spec_hab$Woody.savannas +
      spec_hab$Savannas + spec_hab$Grasslands
  ) / 5
spec_hab$water <-
  rowMeans(spec_hab[, c("Water", "Permanent.wetlands", "Snow.and.ice")], na.rm =
             TRUE)

spec_dist <-
  read.csv(
    "C:\\Users\\xxxx\\Dropbox\\legion2\\host_density_d_distance.csv",
    stringsAsFactors = FALSE
  )
spec_desp <-
  read.csv(
    "C:\\Users\\xxxx\\Dropbox\\legion2\\dispersal_distance_all.csv",
    stringsAsFactors = FALSE
  )
spec_desp$binom <- gsub("_", " ", spec_desp$binom)

###unique
dis1 <- sort(unique(diseases$disease))

x = 0
## change ZIKA TO HUMAN->VECTOR->HUMAN <---- in fact need to relax that generally
for (x in (x + 1):length(dis1)) {
  print(x)
  
  ##choose one disease
  dis2 <- diseases[diseases$disease == dis1[x], ]
  if (length(dis2$tag[!is.na(dis2$tag)]) == 0) {
    print("PROBLEM1")
  }#;next}
  
  ###remove duplicates
  dis2 <- dis2[!duplicated(dis2$name1), ]
  
 ### what is this
  if (paste("", dis1[x], " ALL4.r", sep = "") %in% list.files(
    "C:/Users/Public/Documents/per_disease3/",
    pattern = " ALL4.r",
    full.names = FALSE
  )) {
      next
    }
  
   ##endemic region
  dis_trans <- d1[d1$name == dis1[x], ]
  if (dis_trans$no.countries > 150) {
    next
  }
  ccc <- strsplit(dis_trans$countries, ",")[[1]]
  #if(length(ccc)>=75){print("not widespread");next}
  
  if (is.na(ccc[1])) {
    print("no recognised countries")
    next
  }
  countr <- wrld_simpl[wrld_simpl$ISO2 %in% ccc , ]
  #if(nrow(countr)!=length(ccc)){print("issue with countries");print(ccc)}#;next}
  if (nrow(countr) == 0) {
    print("no recognised countries")
    next
  }
  
  regn <-
    wrld_simpl[wrld_simpl$SUBREGION %in% countr$SUBREGION, ] ##psuedoabsence region
  
  regn2 <-
    wrld_simpl[wrld_simpl$REGION %in% countr$REGION, ] ##psuedoabsence region
  
  ##work out if row is in country
  regn2 <-
    suppressWarnings(crop(regn2, remove_small_islands(regn2, min_val = 10)))
  
  #template2<-fasterize(st_as_sf(regn2),template)
  template2 <- crop(template, regn2)
  
  ### find overlap of scenarios!!
  linkX <- link3[link3$species %in% dis2$name1, ]
  linkX$codeX <- paste(linkX$year, linkX$RCP, sep = "_")
  setDT(linkX)
  setDT(dis2)
  setkey(linkX, species)
  setkey(dis2, name1)
  linkX <- linkX[dis2, allow.cartesian = TRUE]
  linkX <- linkX[!is.na(linkX$codeX), ]
  lv <- unique(linkX$codeX[linkX$type == "vectors"])
  lh <- unique(linkX$codeX[linkX$type == "hosts"])
  if (length(lv[!is.na(lv)]) == 0 | length(lh[!is.na(lh)]) == 0) {
    bothmod <- unique(linkX$codeX)
  } else{
    bothmod <- intersect(lh, lv)
  }
  
  if (length(bothmod) == 0) {
    print("ERROR NO SCENARIO OVERLAP")
    print("PROBLEM2")
    next
  }
  
  ##randomly choose hosts if too many
  if (nrow(dis2[dis2$type == "hosts"]) > 10)	{
    row_drop <-
      sample((1:nrow(dis2))[dis2$type == "hosts"], nrow(dis2[dis2$type == "hosts"]) -
               10, replace = FALSE)
    
    dis2 <- dis2[-row_drop, ]
    
  }
  
  
  ##randomly choose vectors if too many
  if (nrow(dis2[dis2$type == "vectors"]) > 10)	{
    row_drop <-
      sample((1:nrow(dis2))[dis2$type == "vectors"], nrow(dis2[dis2$type == "vectors"]) -
               10, replace = FALSE)
    
    dis2 <- dis2[-row_drop, ]
    
  }
  
  res4 <- NULL
  res5 <- NULL
  seconds_as_hosts = TRUE
  
  ##randomise
  dis2 <- dis2[sample(1:nrow(dis2), nrow(dis2), replace = FALSE), ]
  
  #system.time
  st2 <- Sys.time()

  for (z in 1:nrow(dis2)) {
    #z=z+1
    dis3 <-
      dis2[z, ] #### SORT OUT MULTIPLIE FILES PROBLEM  ### neeed more masks
    
    ### remove GLW dat
    if (dis3$name1 %in% c("ducks",
                          "chickens",
                          "cattle",
                          "ducks",
                          "goats",
                          "human",
                          "pigs",
                          "sheep")) {
      next
    }
    
    ##known locations
    load(file = as.vector(points3[points3$species == dis3$name1, "filen"])) ##called data1$data1
    dt1 <- data.frame(coordinates(data1$data1), data1$data1@data)
    coordinates(dt1) <-  ~ lon + lat
    projection(dt1) <- projection(countr)
    dt1$names2 <- tolower(dt1$name)
    #dt1$names2<-gsub(" ","_",dt1$names2)
    dt1$names2 <- gsub("_", " ", dt1$names2)
    species1 <-
      unique(dt1$names2[!is.na(dt1$names2)])##species in gbif
    ####function to count number of spaces in a string
    no_spaces <-
      function(x, what = " ")
        lengths(regmatches(x, gregexpr(what, x)))
    species1 <- species1[no_spaces(species1) < 2]
    num_spec <- length(species1)
    species1 <-
      species1[sapply(
        species1,
        FUN = function (x)
          is.na(as.numeric(x) == ((as.numeric(
            x
          )) / 1))
      )] ###species in gbif
    if (length(species1) > 0) {
      species2 <-
        unique(read.table(
          text = species1,
          sep = " ",
          stringsAsFactors = FALSE
        )$V1)
    } else{
      species2 <- c()
    } ##genus in gbif
    
    ###habitats etc
    ttt <- gsub("_1_subs_", ";1;", dis3$name1, fixed = TRUE)
    ttt <- gsub("_2_subs_", ";2;", ttt, fixed = TRUE)
    ttt <- gsub("_3_subs_", ";3;", ttt, fixed = TRUE)
    ttt <- gsub("_4_subs_", ";4;", ttt, fixed = TRUE)
    species1a <- strsplit(ttt, ";")[[1]][1]
    species1a <- strsplit(species1a, "_")[[1]][1]
    if (length(strsplit(species1a, " ")[[1]]) > 1) {
      species2a <-
        strsplit(species1a, " ")[[1]][1]
    } else{
      species2a <- species1a
    }
    species2 <- unique(c(species2, species2a))##genus
    species1a <- gsub("_", " ", species1a)
    species1 <- unique(c(species1, species1a))###species
    
    ## add in 
    spec_dist$lower2 <-
      tolower(paste(spec_dist$Taxa, spec_dist$lower, sep = " "))
    
    ## get the species data
    spec_hab2 <- spec_hab[spec_hab$species %in% species1, ]
    spec_dist2 <- spec_dist[spec_dist$lower2 %in% species1, ]
    spec_desp2 <- spec_desp[spec_desp$binom %in% species1, ]
    if (nrow(spec_hab2) == 0) {
      spec_hab2 <- spec_hab[spec_hab$species %in% species2, ]
    }
    if (nrow(spec_dist2) == 0) {
      spec_dist2 <- spec_dist[tolower(spec_dist$Taxa) %in% species2, ]
    }
    if (nrow(spec_desp2) == 0) {
      spec_desp2 <- spec_desp[spec_desp$binom %in% species2, ]
    }
    
    ##put in mean if still none <- ADD MORE TO DATAFRAME
    ##spec_hab ok
    if (nrow(spec_desp2) == 0) {
      spec_desp2 <-
        data.frame(
          binom = NA,
          min = mean(spec_desp$min, trim = 0.1, na.rm = TRUE),
          max = mean(spec_desp$max, trim = 0.1, na.rm = TRUE),
          mean = mean(spec_desp$mean, trim = 0.1, na.rm = TRUE),
          min_mean = mean(spec_desp$min_mean, trim = 0.1, na.rm = TRUE),
          max_mean = mean(spec_desp$max_mean, trim = 0.1, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
    }
    if (nrow(spec_dist2) == 0) {
      spec_dist2 <-
        data.frame(
          Taxa = NA,
          lower = NA,
          Average.of.Mass = mean(
            spec_dist$Average.of.Mass,
            trim = 0.1,
            na.rm = TRUE
          ),
          Average.of.walking.speed = mean(
            spec_dist$Average.of.walking.speed,
            trim = 0.1,
            na.rm = TRUE
          ),
          Average.of.d = mean(
            spec_dist$Average.of.d,
            trim = 0.1,
            na.rm = TRUE
          ),
          Average.of.density = mean(
            spec_dist$Average.of.density,
            trim = 0.1,
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
    }
    
    ###if no mean value make mean value
    if (length(na.omit(spec_desp2$mean)) == 0) {
      spec_desp2$mean = mean(c(spec_desp2$min, spec_desp2$max), na.rm = TRUE)
    }
    if (length(na.omit(spec_desp2$mean)) == 0) {
      print("no dispersal data")
      spec_desp2 <-
        data.frame(
          binom = NA,
          min = mean(spec_desp$min, trim = 0.1, na.rm = TRUE),
          max = mean(spec_desp$max, trim = 0.1, na.rm = TRUE),
          mean = mean(spec_desp$mean, trim = 0.1, na.rm = TRUE),
          min_mean = mean(spec_desp$min_mean, trim = 0.1, na.rm = TRUE),
          max_mean = mean(spec_desp$max_mean, trim = 0.1, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
    }


    #### limit to each species
    link4 <- link3[link3$species == dis3$name1, ]
    #link4$code<-paste(link4$model,link4$year,link4$RCP,sep="_")
    link4$codeX <- paste(link4$year, link4$RCP, sep = "_")
    
    ###remove those with no overlap for hosts and vectors
    ###but what what are best models <-choose both here
    link4 <- link4[link4$codeX %in% bothmod, ]
    link4 <- link4[order(link4$year, link4$RCP), ]
    link4 <-
      link4[c((1:nrow(link4))[link4$RCP == 999], (1:nrow(link4))[!link4$RCP ==
                                                                   999]), ]
    link4$RCP[link4$RCP == 999] <- 2.6
    
    #summarise gas parameters
    link4$speed <-
      mean(spec_dist2$Average.of.walking.speed, na.rm = TRUE)
    link4$d <- mean(spec_dist2$Average.of.d, na.rm = TRUE)
    link4$density <-
      mean(spec_dist2$Average.of.density * num_spec, na.rm = TRUE)
    link4$density[link4$density > 3000000] <- 3000000
    
    
    ##PRESENT DAY FIRST <--------------------------------------------------------START
    link5 <- link4[link4$year == 2010,]
    
    if (nrow(link5) == 0) {
      print(paste0(link4[1, 1], "_maxent.r"))
      print("no present day")
      next
    }
    
    link5 <- link5[sample(1:nrow(link5), 1), ]
    
    ###remove present day
    link4 <- link4[link4$year != 2010, ]
    
    
    ###future climate data
    clims <- raster(link5$filen)

    ##crop
    clims2 <- tryCatch(
      crop(clims, template2),
      error = function(e)
        e
    )
    if ((class(clims2)[1] == "simpleError") == TRUE) {
      t1 <-
        list.files("V:\\new_present2\\", full.names = TRUE)[list.files("V:\\new_present2\\") ==
                                                              paste0(link5$species, "_present_XXX.tif")]
      
      clims <- raster(t1)
      clims2 <- tryCatch(
        crop(clims, template2),
        error = function(e)
          e
      )
    }
    
    if ((class(clims2)[1] == "simpleError") == TRUE) {
      print("Error with climate layer")
    }#;next}
    
    rm(clims)
    clims2 <- mask(clims2, template2)

    ###future land-use data
    futlh <- futlg[futlg$year == link5$year & futlg$RCP == link5$RCP, ]
    lc2 <- stack(futlh$filen)
    names(lc2) <- futlh$variable
    
    ####chose all primary and make not forest
    f1 <- subset(lc2, (1:nlayers(lc2))[futlh$variable == "gothr"])
    nf1 <- f1
    f1[values(fnf) != 1] <- 0
    nf1[values(fnf) == 1] <- 0
    
    ####chose all secondary and make forest/not forest
    s1 <- subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"])
    #wood1<-max(subset(lc2,(1:nlayers(lc2))[futlh$variable %in% c("gfsh1","gfsh2","gfsh3")]),na.rm=TRUE)
    
    ###pasture as grasslands
    ### snow ice water
    water <-
      1 - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable ==
                                                                                                    "gothr"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable ==
                                                                                                                                                                                                          "gcrop"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"])
    #water2<-focal(water)
    
    ## create habitat data for those that at missing it
    if (nrow(spec_hab2) == 0) {
      habsX <- extract(subset(lc2, c(1, 5:8)), dt1, method = 'bilinear')
      meds2 <- apply(habsX, 2, function(x)
        mean(x, trim = 0, na.rm = TRUE))
      meds4 = (meds2 / 2) + 0.5 ### push towards urban due to bias
      #}
      
      lcsuit <-
        (meds4[2] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gothr"])) +
        (meds4[4] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"])) +
        (meds4[3] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"])) +
        (meds4[1] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gcrop"])) +
        (meds4[5] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"]))
      
    } else{
      ###sum up suitability
      #lcsuit<-(wood1*mean(spec_hab2$f,na.rm=TRUE))+(mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
      lcsuit <-
        (f1 * mean(spec_hab2$f, na.rm = TRUE)) + (nf1 * mean(spec_hab2$nf, na.rm =
                                                               TRUE)) + (water * mean(spec_hab2$water, na.rm = TRUE)) + (mean(spec_hab2$Grasslands, na.rm =
                                                                                                                                TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"])) + (mean(spec_hab2$Croplands, na.rm =
                                                                                                                                                                                                            TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gcrop"])) + (mean(spec_hab2$Urban.and.built.up, na.rm =
                                                                                                                                                                                                                                                                                        TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"]))
      meds4 <- NA
    }
    
    
    lcsuit <- resample(lcsuit, clims2)
    
    ##make results table
    res1 <-
      data.table(
        disease = dis3$disease,
        type = dis3$type,
        species = link5$species,
        speed = link5$speed,
        d = link5$d,
        density = link5$density,
        cell.id = values(template2),
        clim = values(clims2),
        lc_suit = values(lcsuit),
        i = 1
      )
    #res1<-res1[!is.na(species),]
    res1 <- res1[!is.na(res1$clim), ]
    res1 <- res1[res1$clim > 0.01, ]
    names(res1)[8:9] <-
      paste(names(res1)[8:9], "present", sep = "_")###lazy sort out!
    res1[, merge1 := paste(cell.id, type, sep = "-")]
    #IP_pres<-IP
    rm(lcsuit, clims2)
    
    ########################################### DO FUTURE	#################################

    ## choose from link4 present day forward
    gc()
    rbl <- function(...)
      rbindlist(list(...))
    cl <- makeCluster(16)
    registerDoParallel(cl)
    
    #st1<-Sys.time()
    res3 = foreach(
      ii = 1:nrow(link4),
      .combine = rbl,
      .packages = c("raster", "data.table")
    ) %dopar%  {

      ##choose one climate model
      link5 <- link4[ii,]
      
      ###future climate data
      clims <- raster(link5$filen)
      clims2 <- tryCatch(
        crop(clims, template2),
        error = function(e)
          e
      )
      if ((class(clims2)[1] == "simpleError") == TRUE) {
        return()
      }
      
      rm(clims)
      clims2 <- mask(clims2, template2)

      ###future land-use data
      futlh <- futlg[futlg$year == link5$year &
                       futlg$RCP == link5$RCP, ]
      lc2 <- stack(futlh$filen)
      names(lc2) <- futlh$variable
      
      ####chose all primary and make not forest
      f1 <- subset(lc2, (1:nlayers(lc2))[futlh$variable == "gothr"])
      nf1 <- f1
      f1[values(fnf) != 1] <- 0
      nf1[values(fnf) == 1] <- 0
      
      ####chose all secondary and make forest/not forest
      s1 <- subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"])
      
      ###pasture as grasslands
      ### snow ice water
      water <-
        1 - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable ==
                                                                                                      "gothr"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable ==
                                                                                                                                                                                                            "gcrop"]) - subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"])
      #wat  ## create habitat data for those that at missing it
      if (nrow(spec_hab2) == 0) {
        ##only create first time around on current land-use
        
        lcsuit <-
          (meds4[2] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gothr"])) +
          (meds4[4] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gsecd"])) +
          (meds4[3] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"])) +
          (meds4[1] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gcrop"])) +
          (meds4[5] * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"]))
        
      } else{
        ###sum up suitability
        #lcsuit<-(wood1*mean(spec_hab2$f,na.rm=TRUE))+(mean(spec_hab2$Grasslands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gpast"]))+(mean(spec_hab2$Croplands,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gcrop"]))+(mean(spec_hab2$Urban.and.built.up,na.rm=TRUE)*subset(lc2,(1:nlayers(lc2))[futlh$variable=="gurbn"]))
        lcsuit <-
          (f1 * mean(spec_hab2$f, na.rm = TRUE)) + (nf1 * mean(spec_hab2$nf, na.rm =
                                                                 TRUE)) + (water * mean(spec_hab2$water, na.rm = TRUE)) + (mean(spec_hab2$Grasslands, na.rm =
                                                                                                                            TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gpast"])) + (mean(spec_hab2$Croplands, na.rm =
                                                                                                                                                                                                              TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gcrop"])) + (mean(spec_hab2$Urban.and.built.up, na.rm =
                                                                                                                                                                                                                                                                                          TRUE) * subset(lc2, (1:nlayers(lc2))[futlh$variable == "gurbn"]))
        
      }
      
      lcsuit <- resample(lcsuit, clims2)
      rm(f1, nf1, water)
      
      ##make results table
      res_foreach <-
        data.table(
          species = link5$species,
          type = dis3$type,
          RCP = link5$RCP,
          year = link5$year,
          type = dis3$type,
          cell.id = values(template2),
          clim_future = values(clims2),
          lc_suit_future = values(lcsuit),
          i = ii
        )
      res_foreach <- res_foreach[!is.na(res_foreach$clim), ]
      res_foreach <- res_foreach[res_foreach$clim > 0.01, ]
      res_foreach[, merge1 := paste(cell.id, type, sep = "-")]
      res_foreach[, cell_by_year_by_RCP := paste(cell.id, year, RCP, sep =
                                                   "_")]
      
      return(res_foreach[!is.na(species)])

    }## END OF FOREACH
    

    stopCluster(cl)
    
    if (is.null(res4)) {
      res4 <- res1
    } else{
      res4 <- rbindlist(list(res4, res1))
    }
    if (is.null(res5)) {
      res5 <- res3
    } else{
      res5 <- rbindlist(list(res5, res3))
    }
    rm(res1, res3)
    gc()
    
    
  }##end of z loop ## z is species/vectors/hosts
  
  ##no present day to compare and not run any so print("PROBLEM");break x
  if (is.null(res4) | is.null(res5)) {
    print("PROBLEM5")
    next
  }
  
  ## REMOVE EFFECT OF SPECIES
  res4b <-
    res4[, .(
      median(clim_present, na.rm = TRUE),
      median(lc_suit_present, na.rm = TRUE),
      median(speed, na.rm = TRUE),
      median(d, na.rm = TRUE),
      mean(density, na.rm = TRUE)
    ), by = .(type, cell.id, merge1)]
  names(res4b)[4:ncol(res4b)] <-
    c(paste(names(res4)[c(8:9, 4:6)], "mean", sep = "_"))
  rm(res4)
 
  ##future by type cellby RCP
   res5c <-
    res5[, .(
      median(clim_future, na.rm = TRUE),
      sd(clim_future, na.rm = TRUE),
      median(lc_suit_future, na.rm = TRUE),
      sd(lc_suit_future, na.rm = TRUE)
    ), by = .(type, cell_by_year_by_RCP, merge1)]
  names(res5c)[4:ncol(res5c)] <-
    c("clim_future_mean",
      "clim_future_sd",
      "lcsuit_future_mean",
      "lcsuit_future_sd")
  rm(res5)#;gc()
  
  ##single host/vector
  if (length(na.omit(res5c$clim_future_sd)) == 0) {
    res5c$clim_future_sd <- 0
    res5c$lcsuit_future_sd <- 0
  }
  
  ##merge present and future
  setkey(res4b, merge1)
  setkey(res5c, merge1)
  res6x <- res5c[res4b]
  res6x <- res6x[, -"i.type"]
  rm(res5c, res4b)
  gc()
  
  ### go to wide
  res6h <- res6x[type == "hosts", ]
  res6v <- res6x[type == "vectors", ]
  rm(res6x)
  gc()
  
  #### if human vec human but has secondary ##rare??
  if (nrow(res6h) == 0) {
    res6h <- res6v
    
    res6v <-
      data.table(
        cell_by_year_by_RCP = unique(res6h$cell_by_year_by_RCP),
        clim_present_mean = 1,
        lc_suit_present_mean = 1,
        clim_future_mean = 1,
        clim_future_sd = 1,
        lcsuit_future_mean = 1,
        lcsuit_future_sd = 1,
        speed_mean = 1,
        d_mean = 1,
        density_mean = 1
      )
  }
  
  #### if human vec human
  if (nrow(res6h) == 0 &
      dis_trans$Type  %in% c("HUMAN->VECTOR->HUMAN->HUMAN", "HUMAN->VECTOR->HUMAN")) {
    res6h <- res6v
    
    res6v <-
      data.table(
        cell_by_year_by_RCP = unique(res6h$cell_by_year_by_RCP),
        clim_present_mean = 1,
        lc_suit_present_mean = 1,
        clim_future_mean = 1,
        clim_future_sd = 1,
        lcsuit_future_mean = 1,
        lcsuit_future_sd = 1,
        speed_mean = 1,
        d_mean = 1,
        density_mean = 1
      )
  }
  
  #### if cattle/livestock are host
  if (nrow(res6h) == 0 &
      all(
        dis2[dis2$type == "hosts", "species"] %in% c(
          "ducks",
          "chickens",
          "cattle",
          "ducks",
          "goats",
          "human",
          "pigs",
          "sheep"
        )
      )) {
    res6h <- res6v
    
    res6v <-
      data.table(
        cell_by_year_by_RCP = unique(res6h$cell_by_year_by_RCP),
        clim_present_mean = 1,
        lc_suit_present_mean = 1,
        clim_future_mean = 1,
        clim_future_sd = 1,
        lcsuit_future_mean = 1,
        lcsuit_future_sd = 1,
        speed_mean = 1,
        d_mean = 1,
        density_mean = 1
      )
  }
  
  #### if true absence - generate dataframe with same RCP year combinations but all zeros
  if (nrow(res6v) < 2 &
      dis_trans$Type %in% c("HOST->HUMAN", "HOST->HUMAN->HUMAN")) {
    res6v <-
      data.table(
        cell_by_year_by_RCP = unique(res6h$cell_by_year_by_RCP),
        clim_present_mean = 1,
        lc_suit_present_mean = 1,
        clim_future_mean = 1,
        clim_future_sd = 1,
        lcsuit_future_mean = 1,
        lcsuit_future_sd = 1,
        speed_mean = 1,
        d_mean = 1,
        density_mean = 1
      )
  }
  
  ### stop if not
  if (nrow(res6v) == 0) {
    res6v <-
      data.table(
        cell_by_year_by_RCP = unique(res6h$cell_by_year_by_RCP),
        clim_present_mean = 1,
        lc_suit_present_mean = 1,
        clim_future_mean = 1,
        clim_future_sd = 1,
        lcsuit_future_mean = 1,
        lcsuit_future_sd = 1,
        speed_mean = 1,
        d_mean = 1,
        density_mean = 1
      )
    print("problem with vectors/hosts")
  }
  
  #### merge column wise for vectors or dummy vectors
  res6v <-
    res6v[, c(
      "cell_by_year_by_RCP",
      "clim_present_mean",
      "lc_suit_present_mean",
      "clim_future_mean",
      "clim_future_sd",
      "lcsuit_future_mean",
      "lcsuit_future_sd",
      "speed_mean",
      "d_mean",
      "density_mean"
    ), with = FALSE]
  names(res6v)[2:ncol(res6v)] <-
    c(paste(names(res6v)[2:ncol(res6v)], "vector", sep = "_"))
  setkey(res6h, cell_by_year_by_RCP)
  setkey(res6v, cell_by_year_by_RCP)
  res6 <- res6h[res6v, nomatch = 0]
  rm(res6h, res6v)
  gc()
  
  ###make year by cell.id by year key
  res6$cell_by_year_by_RCP <- gsub("_2", ";2", res6$cell_by_year_by_RCP)
  res6$cell_by_year_by_RCP <- gsub("0_", "0;", res6$cell_by_year_by_RCP)
  res6[, c("cell.id.2", "year", "RCP") := tstrsplit(cell_by_year_by_RCP, ";")]
  res6[, cell_by_year := paste(cell.id, year, sep = "_")]
  
  ##merge with livestock
  lt3 <-
    lt2[, c("cell_by_year", names(lt2)[tolower(names(lt2)) %in% dis2$name1]), with =
          FALSE]
  
#### ADD SECONDARY TO HOSTS ONLY WHEN HOSTS>0

  if (ncol(lt3) > 1) {
    lt3[, secondary := round(rowSums(lt3[, 2:ncol(lt3)]), 0)]
    ##subset present
    lt3[, c("cell.id", "time") := tstrsplit(cell_by_year, "_")]
    ltp <- lt3[lt3$time == 2010, ]
    ltp <- ltp[, c("cell.id", "secondary")]
    ltp$cell.id <- as.numeric(ltp$cell.id)
    names(ltp)[2] <- "secondary_present"
    setkey(ltp, cell.id)
    ##future
    lt3 <- lt3[, c("cell_by_year", "secondary")]
    names(lt3)[2] <- "secondary_future"
    setkey(res6, "cell_by_year")
    res6 <- res6[lt3, nomatch = 0] 
    ##match present
    setkey(res6, cell.id)
    res6 <- res6[ltp, nomatch = 0]
  } else{
    res6[, secondary_future := 1]
    res6[, secondary_present := 1]
  } ## should be zero because adding to hosts where host + vectors>1

  ##make 0 1 so that it not lost as not key
  res6$secondary_future[res6$secondary_future < 1] <- 1
  res6$secondary_present[res6$secondary_present < 1] <- 1
  
  #### get rid of identical columns
  res6 <-
    res6[, names(res6)[!names(res6) %in% c(
      "type",
      "ID",
      "cell_by_year_by_RCP",
      "cell_by_year",
      "merge1",
      "time",
      "id",
      "cell.id.2",
      "i.id",
      "i.ID",
      "i.cell.id.2",
      "i.time",
      "i.time",
      "i.type",
      "i.cell.id",
      "i.cell.id.1"
    )], with = FALSE]
  
  
  #save present day
  save(res6,
       file = paste(
         "C:/Users/Public/Documents/per_disease3/",
         dis3$disease,
         "ALL4.r"
       ))
  rm(res6)
  gc()
}##end of x loop ## per disease loop
