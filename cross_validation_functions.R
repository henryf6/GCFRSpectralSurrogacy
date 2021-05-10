## ******************************************************************** ##
## cross_validation_functions.R
##
## author: Henry Frye
## data created: March 14, 2017
##
## ******************************************************************** ##


#Libraries needed for non-traditional regression types

#for penalized functional regressions
library(refund)
#for patrial least squartes regressions
library(pls)
#for spatial autocorrelation
library(nlme)

####General Random Cross Validation Function####
CVmaster <- function(data, stat.method = c('pfr','plsr','glm','lm'), formula, seed = 6, k = 10,
                     response, error.method = c('rmse','mse','mae'), family) {
  #allow for custom seed setting for reproducibility
  set.seed(seed)
  
  #assign an id based on the number folds, k
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  list <- 1:k
  
  #create some empty objects for the fold loop to fill
  prediction <- data.frame()
  testsetCopy <- data.frame()
  fold_error <- 0
  
  #Fold loop
  for (i in 1:k){
    #Dividing data into training and test sets based on fold number
    trainingset <- subset(data, id %in% list[-i])
    testset <- subset(data, id %in% c(i))
    
    #run one of the following regressions based on the formula supplied
    #in the function call
    mymodel <- if (stat.method == 'pfr') {
      with(trainingset, pfr(formula = as.formula(formula)))
    } else if (stat.method == 'plsr') {
      plsr(formula = as.formula(formula), data =  trainingset, validation = "CV", parallel = TRUE)
    } else if (stat.method == 'glm') {
      glm(formula  =as.formula(formula), data = trainingset, family = family)
    } else lm(as.formula(formula), data = trainingset)
    
    #put the prediction values of the test data in a temporary dataframe
    temp <- as.data.frame(predict(mymodel, testset, type = 'response'))
    
    
    #make a result 2 column dataframe of the original values from 
    # the fold and the predicted values (from the model of the training set)
    result <- cbind(subset(testset, select = response), temp[,1])
    #names(result) <- c("Actual", "Predicted")
    
    #calculate the error based on one of the method supplied
    fold_error[i] <- if (error.method == 'rmse') {
      sqrt( mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) )
    }  else if (error.method == 'mse') {
      mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) 
    } else mean( abs((result[,1]-result[,2])), na.rm = TRUE ) 
    
  }
  #return a vector with length of the number of folds, k of the error rates
  return(fold_error)
  
}


####Factor Level Cross Validation Function####
#It should do really any factor

# to test this I assigned some group factors to the cars data set
#data(cars)
#cars$colors <- c(rep('blue', 25), rep('red',25))

#fun1 <- function(x, column, fn) {
#  fn(x[,column])
#}
#fun1(df, "B", max)

CVFactor <- function(data, stat.method = c('pfr','plsr','lm'), formula,
                      response, dependent, error.method = c('rmse','mse','mae')) {
  #Creat vector of unique factors
  dependent_factors <- unique(data[,dependent])
  
  #create some empty objects for the fold loop to fill
  prediction <- data.frame()
  testsetCopy <- data.frame()
  fold_error <- 0
  
  #Fold loop
  for (i in 1:length(dependent_factors)){
    #Dividing data into training and test sets based on fold number
    trainingset <- data[which(data[,dependent] != dependent_factors[i]),]
    testset <- data[which(data[,dependent] == dependent_factors[i]),]
    
    #run one of the following regressions based on the formula supplied
    #in the function call
    mymodel <- if (stat.method == 'pfr') {
      with(trainingset, pfr(formula = as.formula(formula)))
    } else if (stat.method == 'plsr') {
      plsr(formula = as.formula(formula), data =  trainingset, validation = "CV", parallel = TRUE)
    } else lm(as.formula(formula), data = trainingset)
    
    #put the prediction values of the test data in a temporary dataframe
    temp <- as.data.frame(predict(mymodel, testset))
    
    
    #make a result 2 column dataframe of the original values from 
    # the fold and the predicted values (from the model of the training set)
    result <- cbind(subset(testset, select = response), temp[,1])
    #names(result) <- c("Actual", "Predicted")
    
    #calculate the error based on one of the method supplied
    fold_error[i] <- if (error.method == 'rmse') {
      sqrt( mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) )
    }  else if (error.method == 'mse') {
      mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) 
    } else mean( abs((result[,1]-result[,2])), na.rm = TRUE ) 
    
  }
  #return a vector with length of the number of folds, k of the error rates
  return(fold_error)
  
}



#dummy <- data.frame(uid =  paste('a',1:10, sep =""),
#                     latitude = runif(10),
#                     longitude = runif(10),
#                     x = rf(10,2,3),
#                     y =  rchisq(10,2))

#((sigma^2)* chol(exp(-as.matrix(dist(geo.dummy))/4)) ) %*% rnorm(10)


####Spatially Dependent Cross Validation Function####

CV_spatial <- function(data, formula, stat.method = c('pfr','plsr','lm'), lat , long , identifier, k, seed1 = 5,
                       response, error.method = c('mse','rmse','mae')) {
  #packages you need: dplyr, nlme, refun, pls
  #STEP 1: Create linear model
  mymodel <- if (stat.method == 'pfr') {
    with(data, pfr(formula = as.formula(formula)))
  } else if (stat.method == 'plsr') {
    plsr(formula = as.formula(formula), data =  data, validation = "CV", parallel = TRUE)
  } else lm(as.formula(formula), data = data)
  
  
  #STEP 2: Create new dataframe with residuals
  data1 <- data.frame(residuals = mymodel$residuals,
                      latitude = data[[lat]],
                      longitude = data[[long]])
  
  #STEP 3: Model Spatial Variance-Covariance matrix using previous model residuals 
  # note: we use a decaying exponential function, (later iterations of function may not have it)
  mymodel.gls <- gls(residuals ~ 1, data=data1, correlation=corExp(form= ~latitude + longitude), method = "REML")
  
  #STEP 4: Create Spatial Distance Matrix (Euclidean Distances only option for now)
  # create empty matrix of n x n dimensions
  spatial   <- matrix(nrow = length(data[[lat]]), ncol = length(data[[long]]),
                      dimnames = list(data[[identifier]] , data[[identifier]]))
  
  # filling matrix with eucliean distances between each point
  #diagonals will be 0, i.e.  distance between the same point is 0
  for(i in 1:length(data[[lat]])) {
    for(j in 1:length(data[[long]])) {
      spatial[i,j] <- sqrt( (data[[lat]][i] - data[[lat]][j])^2 + (data[[long]][i] - data[[long]][j])^2)
    }
  }
  
  #STEP 5 Calculate Spatial Covariance Matrix
  r <- coef(mymodel.gls$modelStruct$corStruct, unconstrained = F)#mymodel.gls$modelStruct$corStruct
  sigma <- mymodel.gls$sigma
  cor_i <- (exp((-spatial)/as.numeric(r)))
  #cov_i <- cor_i*sigma^2
  #all.equal(cov2cor(cov_i), cor_i)
  
  # STEP 6 Multiply Choleski Decomposed Correlation Matrix by random gaussian errors
  set.seed(seed1)
  errs <- rnorm(n= length(data[[lat]]))
  depvec <- chol(cor_i) %*% errs
  
  #STEP 7 Merge Dependent strucutre vector to original data and sorts into groups by the number of folds.
  data$depvec <- depvec
  data$k_enth_tile <- ntile(data$depvec, k)
  list <- 1:k
  
  #STEP 8 create some empty objects for the fold loop to fill
  prediction <- data.frame()
  testsetCopy <- data.frame()
  fold_error <- 0
  
  #STEP 9 Run Folds in a loop
  for (i in 1:k){
    #Dividing data into training and test sets based on fold number
    trainingset <- subset(data, k_enth_tile %in% list[-i])
    testset <- subset(data, k_enth_tile %in% c(i))
    
    #run one of the following regressions based on the formula supplied
    #in the function call
    mymodel <- if (stat.method == 'pfr') {
      with(trainingset, pfr(formula = as.formula(formula)))
    } else if (stat.method == 'plsr') {
      plsr(formula = as.formula(formula), data =  trainingset, validation = "CV", parallel = TRUE)
    } else lm(as.formula(formula), data = trainingset)
    
    #put the prediction values of the test data in a temporary dataframe
    temp <- as.data.frame(predict(mymodel, testset))
    
    
    #make a result 2 column dataframe of the original values from 
    # the fold and the predicted values (from the model of the training set)
    result <- cbind(subset(testset, select = response), temp[,1])
    #names(result) <- c("Actual", "Predicted")
    
    #calculate the error based on one of the method supplied
    fold_error[i] <- if (error.method == 'rmse') {
      sqrt( mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) )
    }  else if (error.method == 'mse') {
      mean( ((result[,1]-result[,2])^2), na.rm = TRUE ) 
    } else mean( abs((result[,1]-result[,2])), na.rm = TRUE ) 
    
  }
  #return a vector with length of the number of folds, k of the error rates
  return(fold_error)
  
}

