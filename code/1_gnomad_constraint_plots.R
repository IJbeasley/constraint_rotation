
# Initial exploration of gnomad constraint data
# in particular, the distribution of p values for different populations

####### set up #######
{
  
library(ggplot2)
library(dplyr)

# Load data (see download_gnomad_data.sh for details)
# using gnomad v2.1.1 
gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# What columns are in the data?
print(names(gnomad_lof))

# Look at the first few rows  and columns of the data
head(gnomad_lof[1:5, 1:5])
}




##################### Part 1.  Investigate the distribution of p values for different populations #####################

# p: The estimated proportion of haplotypes with a pLoF variant. 
# Defined as: 1 -sqrt(no_lofs / defined)


######### summaries of p distributions #########

{

message("\n All populations p distribution")
print(summary(gnomad_lof$p))

message("\n AFR p distribution")
print(summary(gnomad_lof$p_afr))

message("\n NFE p distribution")
print(summary(gnomad_lof$p_nfe))

message("\n FIN p distribution")
print(summary(gnomad_lof$p_fin))

message("\n AMR p distribution")
print(summary(gnomad_lof$p_amr))

}



########### Histograms of p ###########



gnomad_lof  %>%
  ggplot(aes(x = p)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "Distribution of p (for all gnomad v2.1.1)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof   %>%
  ggplot(aes(x = p_afr)) +    
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of p (for AFR)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof   %>%
 ggplot(aes(x = p_nfe)) +    
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of p (for NFE)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof   %>%
 ggplot(aes(x = p_fin)) +    
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of p (for FIN)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof   %>%  
 ggplot(aes(x= p_amr)) +    
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of p (for AMR)",
       x = "p (proportion of haplotypes with a pLoF variant)")

    

####### Histograms of p, for more constrainted genes for different populations ########

# filtering at the median ... 
p_filter = 0.001

gnomad_lof  %>% 
  dplyr::filter(p < p_filter) %>%
  ggplot(aes(x = p)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "All populations (p < 0.001)",
       x = "p (proportion of haplotypes with a pLoF variant)")


gnomad_lof  %>% 
  dplyr::filter(p_afr < p_filter) %>%
  ggplot(aes(x = p_afr)) + 
  geom_histogram(bins = 50) + 
  theme_bw() +
  labs(title = "AFR (p_afr < 0.001)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof  %>% 
  dplyr::filter(p_afr < p_filter) %>% 
  dplyr::pull(p_afr) %>% 
  unique() %>% 
  length()


gnomad_lof  %>% 
  dplyr::filter(p_nfe < p_filter) %>%
  ggplot(aes(x = p_nfe)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "NFE (p_nfe < 0.001)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof %>% 
  dplyr::filter(p_fin < p_filter) %>%
  ggplot(aes(x = p_fin)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "FIN (p_fin < 0.001)",
       x = "p (proportion of haplotypes with a pLoF variant)")

gnomad_lof %>%
 dplyr::filter(p_amr < p_filter) %>%
  ggplot(aes(x = p_amr)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "AMR (p_amr < 0.001)",
       x = "p (proportion of haplotypes with a pLoF variant)")


# Number of unique value for AFR and NFE filtered values 


{
  
p_afr_filtered = gnomad_lof %>% 
                 dplyr::filter(p_afr < p_filter) %>%
                 dplyr::pull(p_afr) 
  

  message("Number of genes with p_afr < ", p_filter)
  
  print(length(p_nfe_filtered))
  
  message("Number of unique p_nfe values < ", p_filter)
  
  print(length(unique(p_nfe_filtered)))
  
  
p_nfe_filtered = gnomad_lof %>% 
                 dplyr::filter(p_nfe < p_filter) %>%
                 dplyr::pull(p_nfe) 

message("Number of genes with p_nfe < ", p_filter)

print(length(p_nfe_filtered))

message("Number of unique p_nfe values < ", p_filter)

print(length(unique(p_nfe_filtered)))

}

# ? Perhaps look at p filtered by overall population constraint, and then look at the distribution of p for each population
# e.g. filter(p < 0.001), then look at p_afr, p_nfe, p_fin, p_amr


################## Plots of delta p (i.e. p_afr - p_nfe, p_fin - p_nfe, and p_amr - p_nfe) ####################

gnomad_lof  = gnomad_lof %>% 
                         dplyr::mutate(delta_p_afr_nfe = p_afr - p_nfe,
                                               delta_p_fin_nfe = p_fin - p_nfe,
                                               delta_p_amr_nfe = p_amr - p_nfe)

{
print("Distribution of delta p AFR NFE")
summary(gnomad_lof$delta_p_afr_nfe)

print("Distribution of delta p FIN NFE")
summary(gnomad_lof$delta_p_fin_nfe)

print("Distribution of delta p AMR NFE")
summary(gnomad_lof$delta_p_amr_nfe)

}

ggplot(gnomad_lof, aes(x = delta_p_afr_nfe)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "Distribution of Delta AFR NFE",
                     x = "Delta p (p_afr - p_nfe)")

ggplot(gnomad_lof, aes(x = delta_p_fin_nfe)) +
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "Distribution of Delta FIN NFE",
                     x = "Delta p (p_fin - p_nfe)")

ggplot(gnomad_lof, aes(x = delta_p_amr_nfe)) +
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "Distribution of Delta AMR NFE",
                     x = "Delta p (p_amr - p_nfe)")

# ? Perhaps look at delta p filtered by overall population constraint, and then look at the distribution of delta p for each population

# ? Perhaps look at delta p - this time measured against the overall population constraint
# i.e.p_nfe - p, p_fin - p, p_afr - p, p_amr - p









############## Part 2:  Investigate the distribution of classic_caf across populations #################

# classic_caf: Sum of allele frequencies of pLoFs in the transcript

########### Histograms of classic_caf ###########



gnomad_lof  %>%
  ggplot(aes(x = classic_caf)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "All populations (classic_caf)",
       x = "classic_caf (Sum of pLoFs allele frequencies)")

gnomad_lof %>% 
  ggplot(aes(x = classic_caf_afr)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "AFR (classic_caf)",
       x = "classic_caf (Sum of pLoFs allele frequencies)")

gnomad_lof %>%
  ggplot(aes(x = classic_caf_nfe)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "NFE (classic_caf)",
         x = "classic_caf (Sum of pLoFs allele frequencies)")

gnomad_lof %>%
  ggplot(aes(x = classic_caf_fin)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "FIN (classic_caf)",
       x = "classic_caf (Sum of pLoFs allele frequencies)")

gnomad_lof  %>% 
  ggplot(aes(x = classic_caf_amr)) + 
  geom_histogram(bins = 100) + 
  theme_bw() +
  labs(title = "AMR (classic_caf)",
       x = "classic_caf (Sum of pLoFs allele frequencies)")




########## Scatter plots of classic_caf (one population classic_caf vs classic_caf_nfe) ##########

 gnomad_lof %>% 
   ggplot(aes(x = classic_caf_nfe, y = classic_caf_afr)) + 
   geom_point() + 
   theme_bw() +
   labs(title = "AFR vs NFE (classic_caf)") 

gnomad_lof %>% 
  ggplot(aes(x = classic_caf_nfe, y = classic_caf_fin)) + 
  geom_point() + 
  theme_bw() +
  labs(title = "FIN vs NFE (classic_caf)") 

gnomad_lof %>%
  ggplot(aes(x = classic_caf_nfe, y = classic_caf_amr)) + 
  geom_point() + 
  theme_bw() +
  labs(title = "AMR vs NFE (classic_caf)")

{
  print("Correlation between classic_caf AFR and NFE")
  cor(gnomad_lof$classic_caf_afr,  gnomad_lof$classic_caf_nfe, method = "spearman",  use =  "complete.obs")

  print("Correlation between classic_caf FIN and NFE")
  cor(gnomad_lof$classic_caf_fin,  gnomad_lof$classic_caf_nfe, method = "spearman",  use =  "complete.obs")

  print("Correlation between classic_caf AMR and NFE")
  cor(gnomad_lof$classic_caf_amr,  gnomad_lof$classic_caf_nfe, method = "spearman",  use =  "complete.obs")
}

##################### scatter plot of p ( different population p vs p_nfe) #####################

{

print("Correlation between p AFR and NFE")
cor(gnomad_lof$p_afr, gnomad_lof$p_nfe, method = "spearman",  use =  "complete.obs")

print("Correlation between p FIN and NFE")
cor(gnomad_lof$p_fin, gnomad_lof$p_nfe, method = "spearman",  use = "complete.obs")

print("Correlation between p AMR and NFE")
cor(gnomad_lof$p_amr, gnomad_lof$p_nfe, method = "spearman",  use = "complete.obs")

}
