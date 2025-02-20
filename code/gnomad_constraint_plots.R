

gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

names(gnomad_lof)


hist(gnomad_lof$classic_caf_afr)



{
  
  afr_no_lof = seq(from = 0, to = 1, by = 0.02) #runif(n = 100000)
  nfe_no_lof =  seq(from = 0, to = 1, by = 0.02)#runif(n = 100000)
  
  diff_df = expand.grid(afr_no_lof, nfe_no_lof)
  
  names(diff_df) = c("afr_no_lof", "nfe_no_lof")
  
  diff_df = diff_df %>% 
    mutate(p_afr = 1 - sqrt(afr_no_lof),
           p_nfe = 1 - sqrt(nfe_no_lof))
  
  # p_afr = 1 - sqrt(afr_no_lof)
  # p_nfe = 1 - sqrt(nfe_no_lof)
  
  diff_df =  diff_df %>% 
    # data.frame(afr_no_lof,
    #                    nfe_no_lof,
    #                    p_afr,
    #                    p_nfe) %>% 
    dplyr::mutate(no_lof_diff = abs(afr_no_lof - nfe_no_lof),
                  p_diff = p_afr - p_nfe)
  
  # Filter points where |x - y| > 0.2
  filtered_df <- diff_df %>% dplyr::filter(no_lof_diff > 0.1)
  
  # Split into two groups
  above_df <- filtered_df %>% dplyr::filter(p_afr > p_nfe)  # Above y = x
  below_df <- filtered_df %>% dplyr::filter(p_afr < p_nfe) # Below y = x
  
  above_df_line = above_df %>% 
    dplyr::group_by(p_nfe) %>% 
    slice_min(p_afr)
  
  below_df_line = below_df %>% 
    dplyr::group_by(p_nfe) %>% 
    slice_max(p_afr)
  
  poly2_above <- lm(p_afr ~ poly(p_nfe, 2, raw = TRUE), data= above_df_line)
  poly2_below <- lm(p_afr ~ poly(p_nfe, 2, raw = TRUE), data= below_df_line)
  
  coef_above <- coef(poly2_above)
  coef_below <- coef(poly2_below)
  
  # Print formulas
  cat("Above y=x: p_afr =", coef_above[1], "+", coef_above[2], "* p_nfe +", coef_above[3], "* p_nfe^2\n")
  cat("Below y=x: p_afr =", coef_below[1], "+", coef_below[2], "* p_nfe +", coef_below[3], "* p_nfe^2\n")
  above_df_line$predicted <- predict(poly2_above, newdata = above_df_line)
  below_df_line$predicted <- predict(poly2_below, newdata = below_df_line)
  
}  


# gnomad_lof %>% 
#   ggplot(aes(x = classic_caf_nfe, y = classic_caf_afr)) + 
#   geom_abline(slope =1, intercept = 0, color = "firebrick")+  
#   geom_point() + 
#   theme_bw() +
#   labs(title = "AFR vs NFE (classic_caf)") 

gnomad_lof %>% 
  ggplot(aes(x = classic_caf_nfe, y = classic_caf_fin)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_point() + 
  theme_bw() +
  labs(title = "FIN vs NFE (classic_caf)") 

##################### Histograms of p 
{

hist(gnomad_lof$p, main = "Distribution of p (for all gnomad v2.1.1)")
summary(gnomad_lof$p)


hist(gnomad_lof$p_afr, main = "Distribution of p (for AFR)")
summary(gnomad_lof$p_afr)

hist(gnomad_lof$p_fin, main = "Distribution of p (for FIN)")
summary(gnomad_lof$p_fin)

hist(gnomad_lof$p_nfe, main = "Distribution of p (for NFE)")
summary(gnomad_lof$p_nfe)

}

##################### Plots of p vs p 

library(ggplot2)
library(dplyr)

gnomad_lof$above_formula <- predict(poly2_above, newdata = gnomad_lof)
gnomad_lof$below_formula <- predict(poly2_below, newdata = gnomad_lof)

filtered_gnomad_lof = gnomad_lof %>% 
  dplyr::filter(p_afr >= above_formula | p_afr <= below_formula)

filtered_gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_afr)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "AFR vs NFE (p)") 

gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_fin)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "FIN vs NFE (p)")


gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_amr)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_point() + 
  theme_bw() +
  labs(title = "AMR vs NFE (p)")

{
cor(gnomad_lof$p_afr, gnomad_lof$p_nfe, method = "spearman",  use =  "complete.obs")

cor(gnomad_lof$p_fin, gnomad_lof$p_nfe, method = "spearman",  use = "complete.obs")

cor(gnomad_lof$p_amr, gnomad_lof$p_nfe, method = "spearman",  use = "complete.obs")

}

################## Plots of delta p

gnomad_lof  = gnomad_lof %>% 
  dplyr::mutate(delta_p_afr_nfe = p_afr - p_nfe,
                delta_p_fin_nfe = p_fin - p_nfe)

summary(gnomad_lof$delta_p_afr_nfe)

summary(gnomad_lof$delta_p_fin_nfe)


hist(gnomad_lof$delta_p_afr_nfe, main = "Distribution of Delta AFR NFE", xlim = c(-0.5, 0.5))

hist(gnomad_lof$delta_p_fin_nfe, main = "Distribution of Delta FIN NFE", xlim = c(-0.5, 0.5))

##########################################