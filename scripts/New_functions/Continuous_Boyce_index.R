#########
# Function to calculate the Boyce index
# Simplified by: https://rdrr.io/cran/modEvA/ 
# Barbosa A.M. [aut], Brown J.A. [aut], Jimenez-Valverde A. [aut], Real R.z [aut], A. Marcia Barbosa [cre]
######
#
#
 
boyce_index<-function (obs = NULL, pred = NULL, n.bins = NA,plot.b=F, 
                        res = 50,plot = TRUE,na.rm = TRUE, ...)
  {

  
  # Extract the predicted values from the raster and the observations
    p.obs <- terra::extract(pred, obs,cells=T) ; p.obs <- p.obs[!duplicated(p.obs$cell),-1]
    p.obs %>% unlist() %>% unname() ; p.obs<-na.omit(p.obs)
    
    fit <- na.omit(terra::values(pred))

  # Functino to extract the parameters of for the boyce calculation
    boycei <- function(interval, p.obs, fit){
                      pi <- sum(as.numeric(p.obs >= interval[1] & p.obs <= interval[2]))/length(p.obs) # Proportion of predicted presences
                      ni <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) # Number of cells
                      ei <- ni/length(fit) # Expected presences if random
                      return(rbind(bin.N = ni, predicted = pi, expected = ei, 
                                   boycei = round(pi/ei, 10) # Proportion of predicted and expected
                      ))
                      }
  
  # Range of values for the intervals and bin widths
  mini <- min(fit, p.obs,na.rm=TRUE)
  maxi <- max(fit, p.obs,na.rm=TRUE)
  
  # Configuration of the intervals and bin widths (values)
  bin.width <- (max(fit) - min(fit))/10 # bin.width
  
  vec.mov <- seq(from = mini, to = maxi - bin.width, 
                 by = (maxi - mini - bin.width)/res)
  
  vec.mov[res + 1] <- vec.mov[res + 1] + 1
  interval <- cbind(vec.mov, vec.mov + bin.width)
    
  # Calculate the Boyce index for the designated moving window
  boycei.result <- t(apply(interval, 1, boycei, p.obs, fit))
  f <- boycei.result[, 4]
  
  to.keep <- which(!is.nan(f))
  f <- f[to.keep]
  
  # Can we calculate the global index
  if (length(f) < 2) {
    b <- NA
  } else {
    r <- 1:length(f)
    j <- vec.mov[to.keep][r]
    b <- cor(f[r], j, method = "spearman")
  }
  
  HS <- apply(interval, 1, sum)/2 # median for each interval
  
  if (length(n.bins) == 1 & is.na(n.bins)) {
    HS[length(HS)] <- HS[length(HS)] - 1
  }
  
  HS <- HS[to.keep]
  
  if(plot.b==TRUE & !is.na(b)){
    plot(HS,f, ylim = c(0, max(f, na.rm = TRUE)),bty="n", 
         xlab = "Prediction class", ylab = "Predicted / expected ratio", cex = 0.5,pch=19,col="black")
  
    lines(HS[r], f[r])
    }

  #
  return(list(bins = data.frame(bin.N = boycei.result[to.keep,1],
                                bin.min = interval[to.keep, 1],
                                bin.max = interval[to.keep,2], 
                                bin.median = HS, 
                                predicted = boycei.result[to.keep,2], 
                                expected = boycei.result[to.keep, 3], 
                                PE.ratio = f),
                                Boyce = b))
      }
