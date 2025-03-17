
# Now let's look at the distribution of p values for the different populations

# And
# Quantify comparisons in gene features for 
# genes that have different p between AFR and NFE groups
# vs genes that do not

{
# Load the data
# using gnomad v2.1.1
gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

library(ggplot2)
library(dplyr)
library(plotly)

}

####### Part i. Get a feature of testis biased expression ##########
{
# Load GTEx data
# Median TPM values for each gene in each tissue
gtex_med_tpm = data.table::fread("data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct")

# Perform PCA on the log transformed data
pca = prcomp(
               gtex_med_tpm %>% 
               select(-Description) %>% 
               tibble::column_to_rownames("Name") %>% 
               as.matrix() %>% 
               +1 %>% 
               log()
             )

}

{
data.frame(pve = pca$sdev/sum(pca$sdev)) %>% 
                  ggplot(aes(x=1:length(pca$sdev), y=pve)) + 
                 geom_col() + 
                 theme_bw() + 
                 labs(x= "PCA Components",
                         y = "Proportion of Variance Explained") 
  
# Compute correlation between PC components with Testis expression
testis_expr <- gtex_med_tpm %>%
                       dplyr::pull(Testis) 

pc_scores <- as.data.frame(pca$x)

pc_cor <- cor(pc_scores, 
                     testis_expr, 
                      use = "pairwise.complete.obs", 
                      method = "spearman"
                      )



# Sort by absolute correlation - to find the PCs that are most correlated with testis expression
pc_cor_df <- data.frame(PC = colnames(pc_scores), Correlation = pc_cor)
pc_cor_df <- pc_cor_df %>% arrange(desc(abs(Correlation)))

message("\n The top 5 PCs most correlated with testis expression")
print(pc_cor_df[1:5,])

# Plot PC1 vs PC6 (the top 2 PCs most correlated with testis expression)
plot = pca$rotation %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Tissue") %>%  
  dplyr::mutate(Testis = ifelse(Tissue == "Testis", 1, 0)) %>% 
  dplyr::mutate(Testis = as.factor(Testis)) %>% 
  ggplot(aes(x=PC1, y=PC6, col = Testis, Tissue = Tissue)) + 
  geom_point(size = 3) + 
  theme_bw() + 
  labs(legend_title = "Testis-biased expression", title = "GTEx Median TPM per tissue") +
  xlab("PC1 (32.98%)") + 
  ylab("PC6 (2.86%)") + 
  theme(axis.title = element_text(size = 20),   # X and Y axis titles
        axis.text = element_text(size = 18),    # X and Y axis labels
        legend.text = element_text(size = 12),  # Legend text
        legend.title = element_text(size = 18), # Legend title
        plot.title = element_text(size = 25, face = "bold") )

print(plot)
ggsave("presentation_figs/pca_plot_pc1_v_pc6.png", width = 15, height = 10, units = "cm")

plotly::ggplotly(plot)
}

{
# Clearly, PC6 is a good distinguishing factor for testis expression
# It separates the testis from other tissues


# How many genes would be distinguished by PC6 as being testis genes?
testis_genes = pca$x %>% 
  as.data.frame() %>% 
  dplyr::select(PC6) %>% 
  dplyr::filter(PC6 < -0.5) %>% 
  row.names()  %>%
  unique() 

message("\n Number of genes that are testis genes")
print(length(testis_genes))
  
message("\n Number of genes in the GTEX dataset")
print(
  pca$x %>% 
  as.data.frame()  %>% 
  row.names() %>% 
  unique() %>% 
  length()
)

# Remove the dot and number after it - to match with the gene names in gnomad data
testis_genes <- sub("\\.\\d+$", "", testis_genes)

# Add a column to the gnomad data to indicate if the gene is a testis gene
gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(is_testis_gene = ifelse(gene_id %in% testis_genes, 1, 0)) 

message("\n Number of testis genes in the gnomad data")  
print(
  gnomad_lof %>% 
  group_by(is_testis_gene) %>% 
  summarise(n = n())
)
}

####################### Define what constitutes a significant difference in p (constraint metric) values ######################

#  Let's define a significant difference in the
#  proportion of haplotypes without a pLoF variant between two populations
# as delta

delta = 0.1

{
set_pop_diff_thresh <- function(delta = 0.1){ 
  
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
  
  #assign("above_formula", above_formula, envir = .GlobalEnv)
  assign("coef_below", coef_below, envir = .GlobalEnv)
  assign("coef_above", coef_above, envir = .GlobalEnv)
  
  
  
}

set_pop_diff_thresh()

gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(above_curve = above_formula(p_nfe)) %>% 
  dplyr::mutate(below_curve = below_formula(p_nfe))

gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(afr_nfe_p_diff = as.factor(ifelse(p_afr >= above_curve | p_afr <= below_curve, 1, 0))) %>%
  dplyr::mutate(amr_nfe_p_diff = as.factor(ifelse(p_amr >= above_curve| p_amr <= below_curve, 1, 0))) %>%
  dplyr::mutate(fin_nfe_p_diff = as.factor(ifelse(p_fin >= above_curve | p_fin <= below_curve, 1, 0)))

}

{
  
  
  # Plot the curve of different p values between AFR and NFE
  # on the gnomad data
  
  gnomad_lof %>% 
    ggplot(aes(x = p_nfe, y = p_afr)) + #colour = afr_nfe_p_diff)) + 
    geom_abline(slope = 1, intercept = 0, color = "firebrick")+  
    geom_point(alpha = 0.25, size = 3) + 
    geom_line(aes(x = p_nfe, y = above_curve), color = "blue", linewidth = 2, linetype = "dashed") +
    geom_line(aes(x = p_nfe, y = below_curve), color = "blue", linewidth = 2, linetype = "dashed") +
    
    theme_bw() +
    ylim(0,1) + 
    xlim(0,1) + 
    theme(axis.title = element_text(size = 20),   # X and Y axis titles
          axis.text = element_text(size = 18),    # X and Y axis labels
          legend.text = element_text(size = 12),  # Legend text
          legend.title = element_text(size = 18), # Legend title
          plot.title = element_text(size = 25, face = "bold") ) + 
    labs(title = "AFR vs NFE (p)")  
  
  
  ggsave("presentation_figs/afr_v_nfe_delta_scatter.png", 
         width = 20, 
         height = 20, 
         units = "cm"
  )  
  
  gnomad_lof %>% 
    ggplot(aes(x = p_nfe, y = p_amr)) + #colour = afr_nfe_p_diff)) + 
    geom_point(alpha = 0.25) + 
    geom_line(aes(x = p_nfe, y = above_curve), color = "blue", linewidth = 0.5, linetype = "dashed") +
    geom_line(aes(x = p_nfe, y = below_curve), color = "blue", linewidth = 0.5, linetype = "dashed") +
    geom_abline(slope =1, intercept = 0, color = "firebrick")+  
    theme_bw() +
    ylim(0,1) + 
    xlim(0,1) + 
    labs(title = "AMR vs NFE (p)") 
  
  gnomad_lof %>% 
    ggplot(aes(x = p_nfe, y = p_fin)) + #colour = afr_nfe_p_diff)) + 
    geom_point(alpha = 0.25) + 
    geom_line(aes(x = p_nfe, y = above_curve), color = "blue", linewidth = 0.5, linetype = "dashed") +
    geom_line(aes(x = p_nfe, y = below_curve), color = "blue", linewidth = 0.5, linetype = "dashed") +
    geom_abline(slope =1, intercept = 0, color = "firebrick")+  
    theme_bw() +
    ylim(0,1) + 
    xlim(0,1) + 
    labs(title = "FIN vs NFE (p)")  
  
}





gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_afr)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = above_line), color = "blue", linewidth = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = below_line), color = "blue", linewidth = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "AFR vs NFE (p)")

  # Now filter for different p values between AFR and NFE
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

# same but for AMR vs NFE
# gnomad_lof %>% 
#   ggplot(aes(x = p_nfe, y = p_amr)) + 
#   geom_abline(slope =1, intercept = 0, color = "firebrick")+  
#   geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
#   geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
#   geom_point(alpha = 0.25) + 
#   theme_bw() +
#   labs(title = "AMR vs NFE (p)")


# gnomad_lof %>% 
#   ggplot(aes(x = p_nfe, y = p_fin)) + 
#   geom_abline(slope =1, intercept = 0, color = "firebrick")+  
#   geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
#   geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
#   geom_point(alpha = 0.25) + 
#   theme_bw() +
#   labs(title = "FIN vs NFE (p)")



############################### Compare genes with different p values between populations ###############################

# Do they have different gene features to those that do not differ between AFR and NFE?

# group genes into different or not between AFR and NFE
gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(diff = ifelse(p_afr >= above_formula | p_afr <= below_formula, 1, 0)) %>% 
  dplyr::mutate(diff = as.factor(diff))

{
message("\n Number of genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(n = n())
)

message("\n Distribution of p_afr and p_nfe for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>%
  summarise(
                      avg_p_afr = mean(p_afr),
                      median_p_afr = median(p_afr),
                      avg_p_nfe = mean(p_nfe),
                      median_p_nfe = median(p_nfe)
                      )
)  
message("\n Distribution of p for not included population,  for genes that differ between AFR and NFE vs those that do not") 
print(
  gnomad_lof %>% 
  group_by(diff) %>%
  summarise(
                      avg_p_amr = mean(p_amr),
                      median_p_amr = median(p_amr),
                      avg_p_fin = mean(p_fin),
                      median_p_fin = median(p_fin),
                      avg_p_oth = mean(p_oth),
                      median_p_oth = median(p_oth),
                      avg_p_sas = mean(p_sas),
                      median_p_sas = median(p_sas),
                      avg_p_eas = mean(p_eas),
                      median_p_eas = median(p_eas)
                      )

)  
message("\n Distrribution of p_afr and p_nfe values, for genes that differ vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(
                     avg_p_afr = mean(p_afr),
                     median_p_afr = median(p_afr),
                     avg_p_nfe = mean(p_nfe),
                     median_p_nfe = median(p_nfe)
                    )


) 

} 

{
message("\n Proportion of testis genes, for those that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(n = n(), 
                    percent_testis_gene = mean(is_testis_gene)
                   )
)
# just below 0.05 significance ... 
print(
  pbinom(q = floor(0.394 * 94), 
       size = 94, 
       prob = 0.328, 
       lower.tail = F
       ) 
)
      
print("\n Distribution of max MAF for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(avg_maf_af = mean(max_af),
                    median_maf_af = median(max_af)
  )
)

message("\n Distribution of gene length features, for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(avg_cds_len = mean(cds_length),
            median_cds_len = median(cds_length),
            avg_gene_len = mean(gene_length),
            median_gene_len = median(gene_length),
            avg_n_exons = mean(num_coding_exons),
            median_n_exons = median(num_coding_exons)
            )
)

message("\n Distribution of Number of LOF sites  for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(
            avg_n_sites = mean(n_sites),
            median_n_sites = median(n_sites)
            )
)
  
} 

########################## subcategories ############# 
# Next step to add to script:
# Create further subcategories of differences 

# 4 categories:
# 1. p_afr > p_nfe, and p_afr and p_nfe are big
# 2. p_afr > p_nfe, and p_afr and p_nfe are small

# 3. p_afr < p_nfe, and p_afr and p_nfe are big
# 4. p_afr < p_nfe, and p_afr and p_nfe are small

# Define big threshold as median
big = 0.001

gnomad_lof  = gnomad_lof  %>% 
                        dplyr::mutate(
                                              unconstrained_p_afr =ifelse(p_afr >= big, 1, 0),
                                              unconstrained_p_nfe = ifelse(p_nfe >= big, 1, 0),
                                              
                                              ) %>% 
                      dplyr:: group_by(diff, unconstrained_p_afr, unconstrained_p_nfe)                        

print("\n Number of genes in each category")
gnomad_lof  %>% 
# group_by(diff, unconstrained_p_afr, unconstrained_p_nfe) %>%
  summarise(n = n())


