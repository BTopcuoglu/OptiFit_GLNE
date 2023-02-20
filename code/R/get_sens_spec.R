library(tidyverse)


infiles <- unlist(snakemake@input)

# infiles <- c(paste0("results/ml/opticlust_denovo/prediction_split_",seq(1,100),".csv"),
#              paste0("results/ml/optifit_self/prediction_split_",seq(1,100),".csv"),
#              paste0("results/ml/optifit_gg/prediction_split_",seq(1,100),".csv"),
#              paste0("results/ml/vsearch_denovo/prediction_split_",seq(1,100),".csv"),
#              paste0("results/ml/vsearch_gg/prediction_split_",seq(1,100),".csv"))

### FUNCTIONS ###########################

get_sensitivity <- function(lookup, x){
  if(x >= max(lookup$specificity)){
    tibble(specificity=x,sensitivity=0)
  } else {
    lookup %>%
      filter(specificity - x > 0) %>%
      top_n(sensitivity, n=1) %>%
      summarize(specificity=x,
                sensitivity = unique(sensitivity))
  }
}

pool_sens_spec <- function(file_name,specificities){
  algorithm <- str_replace(file_name,
                           "results/ml/(.*)/prediction_split_(\\d*).csv", "\\1")
  split <- str_replace(file_name,
                       "results/ml/(.*)/prediction_split_(\\d*).csv", "\\2")
  
  prob <- read_csv(file_name,col_types = cols(Group=col_character(),
                                              dx=col_character(),
                                              .default=col_double()))

  prob_obs <- bind_cols(prob_srn = prob$cancer,
                        observed=prob$dx)

  total <- count(prob_obs, observed) %>%
    pivot_wider(names_from="observed", values_from="n")

  lookup <- prob_obs %>%
    arrange(desc(prob_srn)) %>%
    mutate(is_srn = observed == "cancer") %>%
    mutate(tp = cumsum(is_srn),
           fp = cumsum(!is_srn),
           sensitivity = tp / total$cancer,
           fpr = fp / total$normal) %>%
    mutate(specificity = 1- fpr) %>%
    select(sensitivity, specificity, fpr)

  map_dfr(specificities,get_sensitivity,lookup=lookup) %>%
    mutate(algorithm = algorithm,
           split = split)
}


### VARIABLES ###########################

specificities <- seq(0, 1, 0.01)

### MAIN ################################

all_sens_spec <- map_dfr(infiles,pool_sens_spec,specificities)
write_csv(all_sens_spec,snakemake@output[["outfile"]])
