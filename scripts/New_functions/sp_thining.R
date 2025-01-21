# Function to perform spatial thining based on ecuclidean distances
# adapted from;
#
# Aiello-Lammens, M. E., Boria, R. A., Radosavljevic, A. , Vilela, B. and Anderson, R. P. (2015). spThin: an R package for
# spatial thinning of species occurrence records for use in ecological niche models. Ecography, 38: 541-545. URL
# https://onlinelibrary.wiley.com/doi/10.1111/ecog.01132.
#
#
#' @export thin.algorithm
#' @title Random spatial thining based on euclidean distances
#' 
#' @description \code{sp_thining} applies a distance based randomization approach to
#' spatially thinning occurrence data. 
#' 
#' @param x sf object with the occurrence point occurrence data
#' @param dist.t the distance at which records should be trim, the distance must be specified in the desire units (see units.d)
#' @param units.d the units in which the distances should be measure. e.g "m" for meters of "km" for kilometers (see units package for more details)
#' @param reps The number of times to repeat the thinning process. Each repetition would return a different dataset with the thinning requirements
#' @return records.t: A list object of length 'rep' containing the different trimmed records (sf objects)
#' 

sp_thining<-function(x,
                     dist.t,
                     units.d,
                     reps=1){
  
  list.of.packages<-c("dplyr","sf","units")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Create a list to contain the datasets
  records.t <- vector("list", reps)
  
  y <- x 
  d.y <- y %>% pull(geometry) %>% sf::st_distance()
  units(d.y)<-units::as_units(units.d)
  
  # Remove the units from the distance matrix and convert into a matrix
  d.y <-units::drop_units(d.y) %>% as.matrix()
  d.y<-d.y < dist.t
  
  print(paste0("Distance matrix in ",units.d))
  
  diag(d.y) <- FALSE
  
  ## Set any NA values in the dist matrix to FALSE
  d.y[is.na(d.y)] <- FALSE  
  
  ### Calculate the row sums of the DistMat.save object
  ## This returns the number of elements that are less than
  ## the thin.par for each row
  SumVec.save <- rowSums(d.y)
  
  ## Make a vector of TRUE values of length equal to the number
  ## of rows in the DistMat
  df.keep.save <- rep(TRUE, length(SumVec.save))
  
  for (Rep in seq_len(reps)) {
    ## For each iteration in reps, reset the DistMat and
    ## other indicator variables to original values
    DistMat <- d.y
    SumVec <- SumVec.save
    df.keep <- df.keep.save
    
    ## Perform while loop based on two criteria
    ## 1. The minimum distance between two occurences is less than the 
    ##    thinning parameter 
    ## 2. The number of rows in the resulting data set is greater than 1
    while (any(DistMat) && sum(df.keep) > 1) {
      
      ## Identify the row(s) (occurence) that is within the thin.par distance
      ## to the greatest number of other occurrences. 
      ## If there is more than one row, choose one at random to remove
      RemoveRec <- which(SumVec == max(SumVec))
      if (length(RemoveRec) > 1) {
        RemoveRec <- sample(RemoveRec, 1)
      }
      
      ## Assuming the row chosen above is removed, decrease the 
      ## SumVec object by how many other rows are influenced by its removal
      SumVec <- SumVec - DistMat[, RemoveRec]
      
      ## Set the SumVec value for the row to be removed equal to 0
      SumVec[RemoveRec] <- 0L
      
      ## Set the occ to be ignored in the next iteration of the while loop
      DistMat[RemoveRec, ] <- FALSE
      DistMat[,RemoveRec] <- FALSE
      
      ## Note the occurence for removal from the thinned data set
      df.keep[RemoveRec] <- FALSE
    }
    
    ## Make the new, thinned, data set
    rec.df <- y[df.keep,,drop=FALSE]
    records.t[[Rep]] <- rec.df
  }
  
  ## Order the list object of thinned records by most records
  ## to least
  reduced.rec.order <- unlist(lapply(records.t, nrow))
  reduced.rec.order <- order(reduced.rec.order, decreasing = TRUE)
  records.t <- records.t[reduced.rec.order]
  
  return(records.t)
}

