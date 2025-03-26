


{
  ############################ Load the data ########################
  # using gnomad v2.1.1
  gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
  
  library(ggplot2)
  library(dplyr)
  
}


####################### Define what constitutes a significant difference in p (constraint metric) values ######################

#  Let's define a significant difference in the
#  proportion of haplotypes without a pLoF variant between two populations
# as delta

delta = 0.1

# Function to save function that represents threshold differences in p across ancestry
  set_pop_diff_thresh <- function(delta){ 
    
    pop1_no_lof = seq(from = 0, to = 1, by = 0.01)
    pop2_no_lof =  seq(from = 0, to = 1, by = 0.01)
    
    diff_df = expand.grid(pop1_no_lof, pop2_no_lof)
    names(diff_df) = c("pop1_no_lof", "pop2_no_lof")
    
    diff_df = diff_df %>% 
      dplyr::mutate(
        p_pop1 = 1 - sqrt(pop1_no_lof),
        p_pop2 = 1 - sqrt(pop2_no_lof)
      )
    
    diff_df =  diff_df %>% 
      dplyr::mutate(
        p_diff = p_pop1 - p_pop2,
        no_lof_diff = abs(pop1_no_lof - pop2_no_lof)
      )
    
    filtered_df <- diff_df %>% 
      dplyr::filter(no_lof_diff > delta)
    
    above_df <- filtered_df %>% dplyr::filter(p_pop1 > p_pop2)
    below_df <- filtered_df %>% dplyr::filter(p_pop1 < p_pop2)
    
    above_df_line = above_df %>% 
      dplyr::group_by(p_pop2) %>% 
      dplyr::slice_min(p_pop1)
    
    below_df_line = below_df %>% 
      dplyr::group_by(p_pop2) %>% 
      dplyr::slice_max(p_pop1)
    
    poly2_above <- lm(p_pop1 ~ poly(p_pop2, 2, raw = TRUE), data= above_df_line)
    poly2_below <- lm(p_pop1 ~ poly(p_pop2, 2, raw = TRUE), data= below_df_line)
    
    coef_above <- coef(poly2_above)
    coef_below <- coef(poly2_below)
    
    cat("Above y=x: p_afr =", coef_above[1], "+", coef_above[2], "* p_nfe +", coef_above[3], "* p_nfe^2\n")
    cat("Below y=x: p_afr =", coef_below[1], "+", coef_below[2], "* p_nfe +", coef_below[3], "* p_nfe^2\n")
    
    above_formula <<- function(x_axis_pop) {
      return(coef_above[1] + coef_above[2] * x_axis_pop + coef_above[3] * x_axis_pop ^2)
    }
    
    below_formula <<- function(x_axis_pop) {
      return(coef_below[1] + coef_below[2] * x_axis_pop + coef_below[3] * x_axis_pop ^2)
    }
    
    # put formula in global environment
    assign("coef_below", coef_below, envir = .GlobalEnv)
    assign("coef_above", coef_above, envir = .GlobalEnv)
    
    
    
  }
  
##################### Actually calculate p difference ##############

  # calculate p threshold differences  
  set_pop_diff_thresh(delta = delta)
  
  gnomad_lof = gnomad_lof %>% 
    dplyr::mutate(above_curve = above_formula(p_nfe)) %>% 
    dplyr::mutate(below_curve = below_formula(p_nfe))
  
  # definition: is p different across ancestry?
  gnomad_lof = gnomad_lof %>% 
    dplyr::mutate(afr_nfe_p_diff = as.factor(ifelse(p_afr >= above_curve | p_afr <= below_curve, 1, 0))) %>%
    dplyr::mutate(amr_nfe_p_diff = as.factor(ifelse(p_amr >= above_curve| p_amr <= below_curve, 1, 0))) %>%
    dplyr::mutate(fin_nfe_p_diff = as.factor(ifelse(p_fin >= above_curve | p_fin <= below_curve, 1, 0)))

########################### save p differences dataset ################
  
  gnomad_lof_no_na = gnomad_lof %>% 
    dplyr::filter(!is.na(p))
  
  data.table::fwrite(gnomad_lof_no_na,
                     file = "output/p_diff/gnomad_lof_metrics_delta_p_0.1_na_rm.csv")