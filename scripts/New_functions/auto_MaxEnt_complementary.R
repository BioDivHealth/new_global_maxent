#
############~~~~~~~~~~~~~###############
# AutoMaxent - Complementary functions #
############~~~~~~~~~~~~~###############
#
# a. Create polygons from the coordiantes of a bounding box----
poly_from_ext<-function(x,crs_p=NULL){
  x2<-x[c("xmin", "ymin", "xmax", "ymax")]
  x_double <- as.double(x2)
  names(x_double) <- names(x2)
  class(x_double) <- "bbox"
  
  x_sf <- sf::st_as_sfc(x_double)
  sf::st_crs(x_sf) <- sf::st_crs(crs_p)
  
  return(x_sf)
}

# b. Transform rast into matrix and save the paramters of the raster----
rast_to_vect<-function(x){
  #
  list.of.packages<-c("tidyr","terra","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  #
  dim_r <- c(rows=nrow(x),colums=ncol(x))
  crs_r <- crs(x)
  extend_r <- x %>% terra::ext()
  y <- x %>% terra::as.data.frame(cells=TRUE,na.rm=FALSE)
  
  data <- y[y %>% complete.cases(),]
  NA_index <- y[!y %>% complete.cases(),"cell"]
  
  return(list(dim=dim_r,
              crs=crs_r,
              entent=extend_r,
              tab=data,
              index_missin=NA_index))
}

# c. Variable selection----
var_select<-function(x,
                     VIF.threshold=5,
                     cor.threshold=0.7){
  # prepare the data
  dx <- x
  
  # a. Check Variables VIF values
  select.vars <- vif.sec(yp=dx,threshold=VIF.threshold,
                         return_inf = "select")
  
  # b. Check Variable correlation
  dx <- dx[,colnames(dx) %in% select.vars]
  
  c.int <- cor(dx,use="complete.obs") ; c.int[upper.tri(c.int)]<-NA ; diag(c.int)<-NA
  c.var <- which(abs(c.int) > cor.threshold,arr.ind=TRUE)
  
  if(nrow(c.var)==0){
    print("All correlations are lower than 0.7")
    
  }else{
    pairs.cor <- apply(c.var,1,function(p) data.frame(x=colnames(c.int)[p[1]],y=rownames(c.int)[p[2]])) %>% rbindlist()
    select.vars<-names(dx)[!colnames(dx) %in% c(pairs.cor[,1]%>%unlist())]
  }
  
  return(select.vars)
}

# d. VIF calculation----
vif.sec<-function(yp, # A data.frame or matrix containing the variables from which VIF values need to be calculated
                  vars=names(yp), # A list of vars we want to test
                  threshold=5, # VIF value at which to select the variables
                  return_inf = "full"# "full" or "select", select return the combination of variables that meet the threshold
)
{
  
  # 0. load the needed libraries----
  list.of.packages<-c("usdm")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # Run the calculations
  w<-yp[,colnames(yp)%in%vars]
  w[is.na(w)]<-0
  
  VIF_dat<-list()
  Var_excluded<-list()
  yp<-data.frame(Variables=vars)
  
  # Function to calculate the VIF
  for (i in 1:ncol(w)){
    
    if (i==1){
      yp<-usdm::vif(w)
      
      # Clean up the variables (InF values are coerce to 9000 and NaN variables are exluded)
      yp[is.infinite(yp[,2]),2]<-90000
      null.var <- yp[,1][is.nan(yp[,2])]  ; yp <- yp[!is.nan(yp[,2]),]
      
      y<-yp[yp$VIF==max(yp$VIF),c(1:2)]
      
      # print(y)
      
      Var_excluded[[i]]<-c(as.character(y$Variables))
      #VIF_dat[[i]]<-merge(y,yp,by="Variables",all=TRUE)
      VIF_dat[[i]]<-yp
      
    } else {
      
      if(i %in% c(ncol(w)-1,ncol(w))){
        next()
        
      }else{
        
        A<-w[,-which(names(w) %in% Var_excluded[[i-1]])]
        
        if(length(A)<2){
          Var_excluded[[i]]<-NULL
          break()
          
        }else{
          
          yp<-usdm::vif(A)
          
          # Clean up the variables (InF values are coerce to 9000 and NaN variables are exluded)
          yp[is.infinite(yp[,2]),2]<-90000
          null.var <- yp[,1][is.nan(yp[,2])]  ; yp <- yp[!is.nan(yp[,2]),]
          
          y<-yp[yp$VIF==max(yp$VIF),c(1:2)]
          
          #print(y)
          
          Var_excluded[[i]]<-c(Var_excluded[[i-1]],null.var,as.character(y$Variables))
          #print(y)
          
          VIF_dat[[i]]<-merge(VIF_dat[[i-1]],yp,by="Variables",all=TRUE)
        }
      }
    }
  }
  
  mat.vif<-VIF_dat[[length(VIF_dat)]]
  names(mat.vif)[2:ncol(mat.vif)]<-paste("vif",2:ncol(mat.vif),sep="_")
  
  if(return_inf=="full"){
    return(list(matrix=mat.vif,vars=Var_excluded))
  }
  
  if(return_inf=="select"){
    
    T.a<-apply(mat.vif[,-1],2,FUN=function(x) all(x < threshold,na.rm=TRUE))
    T.a<-names(w)[!names(w) %in% Var_excluded[[which(T.a==TRUE)[1]%>%unname()]]]
    
    return(T.a)
  }
}

# e. Calculate all the possible calculations of a vector
combinations <- function(x, # vector to select
                         choose # Number of elements to consider for each permutation
) {
  
  d <- do.call("expand.grid", rep(list(0:1), length(x)))
  d <-d[rowSums(d) == choose,] ; d <- ifelse(d==1,TRUE,FALSE)
  
  apply(d,1,function(w) x[w]) %>% t()
}

# f. Model Evaluation ----
performance_model<-function(mod.x, # model to test
                            dat, # Test data
                            index_vals,
                            threshold_seq # probability treshold for the classifier
){
  
  y.acc<-data.frame(test=NA)
  
  # Filter the data
  index_data<-complete.cases(dat)
  pa_obs <- index_vals[index_data]
  mod.predictions <- predict(mod.x,dat[index_data,],type="response")
  
  for(l in 1:length(threshold_seq)){
    
    # Create confusion matrix
    Conf.Mat<-table(predict=ifelse(mod.predictions >= threshold_seq[l], 1, 0),real=pa_obs)
    
    if(nrow(Conf.Mat)<2){
      next()
    }
    
    Tn<-Conf.Mat[1,1] ; Fn<-Conf.Mat[1,2]
    Fp<-Conf.Mat[2,1] ; Tp<-Conf.Mat[2,2]
    
    
    A<-list(Accuracy=(Tn+Tp)/(Tn+Tp+Fp+Fn),
            TypeI=(Fp)/(Tn+Tp+Fp+Fn),
            TypeII=(Fn)/(Tn+Tp+Fp+Fn),
            TNR=(Tn)/(Tn+Fp),
            TPR=(Tp)/(Tp+Fn),
            kappa=(2*(Tp*Tn-Fn*Fp))/((Tp+Fp)*(Fp+Tn)+(Tp+Fn)*(Fn+Tn)),
            TSS=((Tp)/(Tp+Fn))+((Tn)/(Tn+Fp)-1)
    )
    
    b <- A %>% unlist() %>% as.data.frame() 
    b$test<- names(A)
    
    y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
    
    #rm(A,b)
  }
  
  if(length(y.acc)<2){
    print("No performance information for the Model")
    return(NULL)
  }else{
    y.acc<-y.acc[!is.na(y.acc$test),] ; rownames(y.acc)<-y.acc$test
    y.acc<-y.acc[,-1] %>% as.matrix()
    
    return(y.acc)
  }
}

# e. Calculation of the Boyce index ----
# Function to calculate the Boyce index
# Simplified by: https://rdrr.io/cran/modEvA/ 
# Barbosa A.M. [aut], Brown J.A. [aut], Jimenez-Valverde A. [aut], Real R.z [aut], A. Marcia Barbosa [cre]
#

boyce_index<-function (obs = NULL, pred = NULL, n.bins = NA, 
                       res = 100,plot = TRUE,na.rm = TRUE, ...)
{
  
  # Extract the predicted values from the raster and the observations
  p.obs <- terra::extract(pred, obs) %>% unlist() %>% unname() ; p.obs<-na.omit(p.obs)
  fit <- na.omit(terra::values(pred))
  
  
  boycei <- function(interval, p.obs, fit) {
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
    b <- cor(f[r], vec.mov[to.keep][r], method = "spearman")
  }
  
  HS <- apply(interval, 1, sum)/2 # median for each interval
  
  if (length(n.bins) == 1 & is.na(n.bins)) {
    HS[length(HS)] <- HS[length(HS)] - 1
  }
  
  HS <- HS[to.keep]
  
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

# End of the script
