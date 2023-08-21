

makePreds <- function(model, testData, trainData){
  
  # This function is far from complete.  Plan is to combine multiple feature matrices
  # that are in different files (testData is a character vector of file names) and concatenate
  # them into 1 feature matrix and then make predictions using the model.
  # First step should be to read in each file individually and concatenate them.
  # But what do I merge on? For some sets we have duplicate DTXSIDs and SMILES.
  
  # Read in and combine all test data
  for (i in 1:length(testData)){
    temp <- read.delim(file = testData[i], header = TRUE, sep = "\t")
    if (i != 1){
      feature_flag <- length(intersect(colnames(temp) == colnames(trainData)) > 5)  # Same feature set
    }
  }
  
  # Limit test data to the features in the model
  test <- testData[,colnames(testData) %in% colnames(trainData)]
  test <- apply(test, 2, as.character)
  test <- as.data.frame(test)
  
  # Predict primary elimination route
  test_class <- predict(rf_classifier, newdata = test, type = "prob")
  
  # Examine results
  sum(test_class[,2] > 0.5)
  hist(test_class[,1])
  hist(test_class[,2])
  
  # Save these results
  TSCAactive_predRoute <- testData[,c(1:3,5,7)]
  TSCAactive_predRoute$probDiffusion <- test_class[,1]
  TSCAactive_predRoute$probFiltration <- test_class[,2]
  TSCAactive_predRoute$class <- "diffusion"
  TSCAactive_predRoute$class[TSCAactive_predRoute$probFiltration > 0.5] <- "filtration"
  
}


#' SimilarityThreshold
#'
#' Determine similarity threshold for domain of applicability calculation
#'
#' For each classification, \code{SimilarityThreshold} determines what
#' distance chemical in a predicted set can be from the training set in order
#' for that prediction to be in the model's domain of applicability
#'
#' @param train \code{data.frame} of training descriptors
#' @param sim_cut fraction from 0 to 1 for how similar records should be, passed to SimilarityThreshold
#' @param method A value for the cutoff distance
#'
#' @return
#' @export
#'
#' @examples SimilarityThreshold(train,sim_cut=0.5,method="jaccard")
#'
SimilarityThreshold <- function(train, sim_cut=0.5, method="euclidean"){
  
  ## Check that necessary package is installed
  #if (!("proxy" %in% rownames(installed.packages()))){
  #  print("Installing proxy packages -- required for this function")
  # install.packages("proxy", "/data/home/kisaacs1/anaconda3/lib/R/library/", repos="http://cran.us.r-project.org")
  #}
  
  ## Pair-wise distance matrix for the training set
  d_matrix <- proxy::dist(train,method=method)
  
  ## Calculate the distance matrix of all descriptors in the training set
  d_matrix <- as.matrix(d_matrix)
  
  ## Calculate the number of nearest neighbors that makes up the average and
  ## standard deviation of the nearest neighbors distances
  ## Wishlist -- re-find reference for using the .35% value, and make it a
  ## function parameter
  knn <- ceiling(0.035*dim(d_matrix)[1])
  
  ## The diagonal of the matrix is not needed, so make it null
  diag(d_matrix) <- NA
  
  ## Now get the nearest knn distances between a molecule and all other
  ## molecules in the training set.
  d_knn <- apply(d_matrix,1,function(x){sort(x)[1:knn]})
  
  ## R calculates the sample standard deviation, I want the population standard
  ## deviation
  std <- sqrt(mean(d_knn ^ 2) - mean(d_knn)*mean(d_knn))
  
  ## Now calculate the similarity threshold
  d_cut <- sim_cut * mean(d_knn) + std
  
  return(d_cut)
}


#' DomainOfApplicability
#'
#' Predict if a record is in a model's domain of applicability
#'
#' For each record in a testing set of data, \code{DomainOfApplicability}
#' determines if that record lies within the domain of applicability of the
#' training set for a given class in a classification model.
#'
#' References:
#' 1. Rational selection of training and test sets for the development of
#'    valid QSAR models, A. Golbraikh, et al., J. Comp.-Aided Mol. Design,
#'    17 (2003), 241 -- 253.
#' 2. Predictive QSAR modeling workflow, model applicability domains, and
#'    virtual screening, A. Tropsha and A. Golbraikh, Curr. Pharm. Design,
#'    13 (2007), 3494 -- 3504.
#'
#' Parameters:
#' @param training \code{data.frame} of training data set that was used to build models
#' @param predicted \code{data.frame} of testing data set that was predicted with models
#' @param train_class name of column in \code{training} that was use for model training
#' @param pred_class name of column in \code{predicted} that was was predicted
#' @param model_list list of models, or "classes", for which to determine domain
#' @param sim_cut fraction from 0 to 1 for how similar records should be, passed to SimilarityThreshold
#' @param method distance metric to use for similarity matrix (based on \code{proxy} package)
#'
#' @return A vector of indices of records in domain of applicability
#'
#' @examples DomainOfApplicability(training_set,testing_set,model_list,metric="jaccard")
#'
DomainOfApplicability <- function(training,predicted,sim_cut=0.5,method="euclidean"){
  
  ## Check that necessary package is installed
  #if (!("proxy" %in% rownames(installed.packages()))){
  #  print("Installing proxy packages -- required for this function")
  # install.packages("proxy", "/data/home/kisaacs1/anaconda3/lib/R/library/", repos="http://cran.us.r-project.org")
  #}
  
  ## Check that train_class and pred_class are indeed columns in their
  ## respective data sets
  # if (!(train_class %in% colnames(training))){
  #   stop(paste("Column ''",train_class,"' not a column in training.",sep="",collapse=""))
  # }
  # if (!(pred_class %in% colnames(predicted))){
  #   stop(paste("Column ''",pred_class,"' not a column in predicted.",sep="",collapse=""))
  # }
  
  
  if (!("DTXSID" %in% colnames(training))){
    print(paste("Warning: Column 'DTXSID' is not a column in 'training', if CASRNs are in a\ndifferent column, they will be included as descriptor."))
  }
  if (!("DTXSID" %in% colnames(predicted))){
    print(paste("Warning: Column 'DTXSID' is not a column in 'predicted', if CASRNs are in a\ndifferent column, they will be included as descriptor."))
  }
  
  ## Only keep descriptors of training data
  cols <- which(!(colnames(training) %in% c("DTXSID")))
  train <- training[,cols]
  
  ## Only keep descriptors of predicted data
  cols <- which(!(colnames(predicted) %in% c("DTXSID")))
  preds <- predicted[,cols]
  
  ## This is the similarity threshold that Tropsha uses in most of his QSAR
  ## validation work: d_cut <- Z*d_avg + d_std
  d_cut <- SimilarityThreshold(train=train, method=method)
  
  ## Calculate the pair-wise distances between molecules in the training set
  ## and predictions
  d_train_preds <- proxy::dist(train, preds, method=method)
  
  ## Now that you know the distances see if they meet the threshold requirement
  ## Find the closed training chemical to the predicted chemicals, and store
  ## that distance
  d_ij <- apply(d_train_preds, 2, "min")
  
  ## Now for each predicted chemical, see if that closest distance is less than
  ## the threshold distance
  in_domain <- sapply(d_ij, function(x){x < d_cut})
  
  ## If it is, then store the index, if not then its not in the domain of
  ## applicability and is not a valid prediction
  idx <- which(in_domain==TRUE)
  
  ## Return data frame indices of records within domain of applicability
  return(idx)
}

