library(tidyverse)

data <- read_csv(snakemake@input[["merged_performance"]])
#data <- read_csv("results/ml/summary/merged_performance.csv")

outfile <- snakemake@output[["outfile"]]

### FUNCTIONS

# calculate the difference in the mean of the metric for the two groups
# sub_data: subset of the performance data for just two groups
# group_name: name of column with group variable
# metric: metric to compare
get_difference <- function(sub_data,group_name,metric){
  # get the mean metric value for each group
  means <- sub_data %>%
    group_by(.data[[group_name]]) %>%
    summarise(meanVal = mean(.data[[metric]]),.groups="drop") %>%
    pull(meanVal)
  # find the difference in the mean between the two groups
  abs(diff(means))
}

# shuffle the groups across the data
# sub_data: subset of the performance data for just two groups
# group_name: column name to shuffle
shuffle_group <- function(sub_data,group_name){
  group_vals <- sub_data %>% 
    pull( {{ group_name }} )
  group_vals_shuffled <- base::sample(group_vals)
  
  data_shuffled <- sub_data %>% 
    mutate( !!group_name := group_vals_shuffled)
  
  return(data_shuffled)
}

# data: the concatenated performance (from what output?)
# metric: metric to compare, either AUC or cv_metric_AUC
# group_name: column with group variables to compare
# group_1: name of one group to compare
# group_2: name of other group to compare
perm_p_value <- function(data, metric, group_name, group_1, group_2, nperm = 10000){
  # check that the metric, group, and both group labels exist in data
  assertthat::has_name(data,metric)
  assertthat::has_name(data,group_name)
  ### FIX THIS check that group_1 and group_2 exist in the data
  # is this necessary? 
  #assertthat::has_name(data %>% distinct({{group_name}}),{{ group_1 }})
  
  # subset results to select metric and group columns and
  # filter to only the two groups of interest
  sub_data <- data %>%
    select({{ metric }},{{ group_name }}) %>%
    filter( .data[[group_name]] == {{group_1}}  | .data[[group_name]] == {{group_2}})
  
  # observed difference
  # quantify the absolute value of the difference in metric between the two groups (observed difference)
  obs <- get_difference(sub_data,{{group_name}},{{metric}})
  
  # shuffled difference
  # quantify the absolute value of the difference in metric between the two groups after shuffling group labels
  metric.null <-  replicate(nperm,get_difference(shuffle_group(sub_data,group_name),group_name,metric))
  
  # n = number of shuffled calculations
  n <- length(metric.null)
  # r = replications at least as extreme as observed effect
  r <- sum(abs(metric.null) >= obs)
  
  # compute Monte Carlo p-value with correction (Davison & Hinkley, 1997)
  p.value=(r+1)/(n+1)
  return(p.value)
}

# data: the concatenated performance (from what output?)
# metric: metric to compare, either AUC or cv_metric_AUC
# group_name: column with group variables to compare
compare_models <- function(merged_data,metric,group_name){
  
  # identify all unique groups in group variable
  groups <- merged_data %>% 
    pull( {{group_name}} ) %>% 
    unique()
  
  # create a table with all possible comparisons of groups
  # without repeating pairings
  p.table <- expand.grid(x=1:length(groups),
                         y=1:length(groups)) %>% 
    filter(x < y) %>% 
    mutate(group1 = groups[x],
           group2 = groups[y])  %>% 
    select(-x, -y) %>% 
    group_by(group1,group2) %>% 
    summarize(p_value = perm_p_value(merged_data,metric,group_name,group1,group2),
              .groups = "drop")
  
  return(p.table)

}

### CALCULATE P-VALUES

pvals <- bind_rows(c(metric="cv_metric_AUC",compare_models(data,"cv_metric_AUC","algorithm")),
                   c(metric="AUC",compare_models(data,"AUC","algorithm"))) 

write_csv(pvals,outfile)
