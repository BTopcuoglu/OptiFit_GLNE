
# Author: Begum Topcuoglu
# Date: 2019-01-14
######################################################################
# Description:
# This script trains and tests the model according to proper pipeline
######################################################################

######################################################################
# Dependencies and Outputs:
#    Model to put to function:
#       1. "L2_Logistic_Regression"
#       2. "L2_Linear_SVM"
#       3. "RBF_SVM"
#       4. "Decision_Tree"
#       5. "Random_Forest"
#       6. "XGBoost"
#    Dataset to put to function:
#         Features: Hemoglobin levels and 16S rRNA gene sequences in the stool
#         Labels: - Colorectal lesions of 490 patients.
#                 - Defined as cancer or not.(Cancer here means: SRN)
#
# Usage:
# Call as source when using the function. The function is:
#   pipeline(data, model)

# Output:
#  A results list of:
#     1. cvAUC and testAUC for 1 data-split
#     2. cvAUC for all hyper-parameters during tuning for 1 datasplit
#     3. feature importance info on first 10 features for 1 datasplit
#     4. trained model as a caret object
######################################################################

######################################################################
#------------------------- DEFINE FUNCTION -------------------#
######################################################################


pipeline <- function(dataset, test, model, outcome=NA, hyperparameters=NULL){

  # -----------------------Get outcome variable----------------------------->
  # If no outcome specified, use first column in dataset
  if(is.na(outcome)){
    outcome <- colnames(dataset)[1]
  }else{
    # check to see if outcome is in column names of dataset
    if(!outcome %in% names(dataset)){
      stop(paste('Outcome',outcome,'not in column names of dataset.'))
    }
  }

  # ------------------Pre-process the full Dataset------------------------->
  # We are doing the pre-processing to the full dataset and then splitting 80-20
  # Scale all features between 0-1
  preProcValues <- preProcess(dataset, method = "range")
  dataTransformed <- predict(preProcValues, dataset)

  test <- test %>%
    select(-Group)

  testTransformed <- predict(preProcValues, test)
  # ----------------------------------------------------------------------->

  # Get outcome variables
  first_outcome = as.character(dataset[,outcome][1])
  outcome_vals = unique(dataset[,outcome])
  if(length(outcome_vals) != 2) stop('A binary outcome variable is required.')
  second_outcome = as.character(outcome_vals[!outcome_vals == first_outcome])
  print(paste(c('first outcome:','second outcome:'),c(first_outcome,second_outcome)))


  # -------------Define hyper-parameter and cv settings-------------------->
  # Define hyper-parameter tuning grid and the training method
  # Uses function tuning_grid() in file ('code/learning/model_selection.R')
  tune <- tuning_grid(dataTransformed, model, outcome, hyperparameters)
  grid <- tune[[1]]
  method <- tune[[2]]
  cv <- tune[[3]]
  # ----------------------------------------------------------------------->

  # ---------------------------Train the model ---------------------------->
  # ------------------------------- 1. -------------------------------------
  # - We train on the 80% of the full dataset.
  # - We use the cross-validation and hyper-parameter settings defined above to train
  # ------------------------------- 2. -------------------------------------
  # We use ROC metric for all the models
  # To do that I had to make changes to the caret package functions.
  # The files 'data/caret_models/svmLinear3.R and svmLinear5.R are my functions.
  # I added 1 line to get Decision Values for linear SVMs:
  #
  #           prob = function(modelFit, newdata, submodels = NULL){
  #             predict(modelFit, newdata, decisionValues = TRUE)$decisionValues
  #           },
  #
  # This line gives decision values instead of probabilities and computes ROC in:
  #   1. train function with the cross-validataion
  #   2. final trained model
  # using decision values and saves them in the variable "prob"
  # ------------------------------- 3. --------------------------------------
  # - If the model is logistic regression, we need to add a family=binomial parameter.
  # - If the model is random forest, we need to add a ntree=1000 parameter.
  #         We chose ntree=1000 empirically.
  # ----------------------------------------------------------------------->
  # Make formula based on outcome
  f <- as.formula(paste(outcome, '~ .'))
  print('Machine learning formula:')
  print(f)
  # Start walltime for training model
  tic("train")
  if(model=="L2_Logistic_Regression"){
  print(model)
  trained_model <-  train(f, # label
                          data=dataTransformed, #total data
                          method = method,
                          trControl = cv,
                          metric = "ROC",
                          tuneGrid = grid,
                          family = "binomial")
  }
  else if(model=="Random_Forest"){
      print(model)
      trained_model <-  train(f,
                              data=dataTransformed,
                              method = method,
                              trControl = cv,
                              metric = "ROC",
                              tuneGrid = grid,
                              ntree=1000) # not tuning ntree
  }
  else{
  print("Did not define a model algorithm to use!")
  }

  # ------------- Output the cvAUC ---------------------------------------------------->
  # Mean cv AUC value over repeats of the best cost parameter during training
  cv_auc <- getTrainPerf(trained_model)$TrainROC
  # ---------------------------------------------------------------------------------->

  # -------------------------- Predict the held-out sample---------------------------->

  # Calculate the test-auc for the actual pre-processed held-out data
  rpartProbs <- predict(trained_model, testTransformed, type="prob")

  # ---------------------------------------------------------------------------------->

  # Calculate the best decision threshold for this model
  library(pROC)
  test_roc <- roc(ifelse(testTransformed[,outcome] == first_outcome, 1, 0), rpartProbs[[1]])
  thr <- coords(test_roc, "best", ret = "threshold")

  # ----------------------------Save metrics as vector ------------------------------->
  # Return all the metrics
  results <- list(cv_auc, rpartProbs, trained_model, thr)
  return(results)
}
