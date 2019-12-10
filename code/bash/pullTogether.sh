#!/bin/bash

# Author: Begum Topcuoglu
# Date: 2018-02-13
#

SEARCH_DIR=data/temp
FINAL_DIR=data/process

samples=`awk '{print $1}' data/process/baxter/glne007.files`


for sample in $samples
do
  	head -1 $SEARCH_DIR/cv_results_"$sample".csv  > $FINAL_DIR/combined_cv_results.csv ; tail -n +2 -q $SEARCH_DIR/cv_results_"$sample".csv >> $FINAL_DIR/combined_cv_results.csv

  	head -1 $SEARCH_DIR/prediction_results_"$sample".csv  > $FINAL_DIR/combined_prediction_results.csv ; tail -n +2 -q $SEARCH_DIR/prediction_results_"$sample".csv >> $FINAL_DIR/combined_prediction_results.csv
done
