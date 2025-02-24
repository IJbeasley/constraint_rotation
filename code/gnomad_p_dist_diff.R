
# Now let's look at the distribution of p values for the different populations

# And
# Quantify comparisons in gene features for 
# genes that have different p between AFR and NFE groups
# vs genes that do not

# Load the data
# using gnomad v2.1.1
gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

library(ggplot2)
library(dplyr)
library(plotly)

####### Part i. Get a feature of testis biased expression ##########

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

print("\n The top 5 PCs most correlated with testis expression")
print(pc_cor_df[1:5,])

# Plot PC1 vs PC6 (the top 2 PCs most correlated with testis expression)
plot = pca$rotation %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Tissue") %>%  
  dplyr::mutate(Testis = ifelse(Tissue == "Testis", 1, 0)) %>% 
  dplyr::mutate(Testis = as.factor(Testis)) %>% 
  ggplot(aes(x=PC1, y=PC6, col = Testis, Tissue = Tissue)) + 
  geom_point() + 
  theme_bw() 


plotly::ggplotly(plot)

# Clearly, PC6 is a good distinguishing factor for testis expression
# It separates the testis from other tissues


# How many genes would be distinguished by PC6 as being testis genes?
testis_genes = pca$x %>% 
  as.data.frame() %>% 
  dplyr::select(PC6) %>% 
  dplyr::filter(PC6 < -0.5) %>% 
  row.names()  %>%
  unique() 

print("\n Number of genes that are testis genes")
length(testis_genes) 
  
print("\n Number of genes in the GTEX dataset")
pca$x %>% 
  as.data.frame()  %>% 
  row.names() %>% 
  unique() %>% 
  length()

# Remove the dot and number after it - to match with the gene names in gnomad data
testis_genes <- sub("\\.\\d+$", "", testis_genes)

# Add a column to the gnomad data to indicate if the gene is a testis gene
gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(is_testis_gene = ifelse(gene_id %in% testis_genes, 1, 0)) 

print("\n Number of testis genes in the gnomad data")  
gnomad_lof %>% 
  group_by(is_testis_gene) %>% 
  summarise(n = n())


####################### Define what constitutes a significant difference in p (constraint metric) values ######################

#  Let's define a significant difference in the
#  proportion of haplotypes without a pLoF variant between two populations
# as delta

delta = 0.1

{
# Under this definition, what differences in p (estimated as 1- sqrt(no_lofs/defined)) are significant?
# Let's estimate it by simulating: 

# simulate the proportion of haplotyes without a pLoF variant in AFR and NFE
afr_no_lof = seq(from = 0, to = 1, by = 0.02)
nfe_no_lof =  seq(from = 0, to = 1, by = 0.02)

# assume every combination of afr_no_lof and nfe_no_lof is possible
diff_df = expand.grid(afr_no_lof, nfe_no_lof)
names(diff_df) = c("afr_no_lof", "nfe_no_lof")

# calculate p_afr and p_nfe values
  diff_df = diff_df %>% 
    dplyr::mutate(p_afr = 1 - sqrt(afr_no_lof),
                           p_nfe = 1 - sqrt(nfe_no_lof)
                           )

# calculate the difference in p_afr and p_nfe values
  diff_df =  diff_df %>% 
                  dplyr::mutate(
                                         p_diff = p_afr - p_nfe
                                         no_lof_diff = afr_no_lof - nfe_no_lof
                                         )

  # Filter points where |afr_no_lof  -  nfe_no_lof| > delta
  filtered_df <- diff_df %>% 
                        dplyr::filter(no_lof_diff > delta)

  # Split into two groups
  above_df <- filtered_df %>% dplyr::filter(p_afr > p_nfe)  # Above y = x (p_afr > p_nfe)
  below_df <- filtered_df %>% dplyr::filter(p_afr < p_nfe) # Below y = x (p_afr < p_nfe)

# Calculate the curve that fits the edge of the above and below points
above_df_line = above_df %>% 
    dplyr::group_by(p_nfe) %>% 
    dplyr::slice_min(p_afr)
  
  below_df_line = below_df %>% 
    dplyr::group_by(p_nfe) %>% 
    dplyr::slice_max(p_afr)
  
poly2_above <- lm(p_afr ~ poly(p_nfe, 2, raw = TRUE), data= above_df_line)
poly2_below <- lm(p_afr ~ poly(p_nfe, 2, raw = TRUE), data= below_df_line)

 coef_above <- coef(poly2_above)
coef_below <- coef(poly2_below)
  
# Print formulas
cat("Above y=x: p_afr =", coef_above[1], "+", coef_above[2], "* p_nfe +", coef_above[3], "* p_nfe^2\n")
cat("Below y=x: p_afr =", coef_below[1], "+", coef_below[2], "* p_nfe +", coef_below[3], "* p_nfe^2\n")

# Now get the predicted values for the above and below points
# To plot the curve that fits the edge of the above and below points
gnomad_lof$above_formula <- predict(poly2_above, newdata = gnomad_lof)
gnomad_lof$below_formula <- predict(poly2_below, newdata = gnomad_lof)

}

# Plot the curve of different p values between AFR and NFE
# on the gnomad data
gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_afr)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "AFR vs NFE (p)")

  # Now filter for different p values between AFR and NFE
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
gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_amr)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "AMR vs NFE (p)")


gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_fin)) + 
  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
  geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.25) + 
  theme_bw() +
  labs(title = "FIN vs NFE (p)")



############################### Compare genes with different p values between populations ###############################

# Do they have different gene features to those that do not differ between AFR and NFE?

# group genes into different or not between AFR and NFE
gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(diff = ifelse(p_afr >= above_formula | p_afr <= below_formula, 1, 0)) %>% 
  dplyr::mutate(diff = as.factor(diff))

print("\n Number of genes that differ between AFR and NFE")
gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(n = n())

print("\n Number of testis genes that differ between AFR and NFE")
gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(n = n(), 
                    percent_testis_gene = mean(is_testis_gene)
                   )

# just below 0.05 significance ... 
pbinom(q = floor(0.394 * 94), 
       size = 94, 
       prob = 0.328, 
       lower.tail = F
       ) 

print("\n Distribution of features for genes that differ between AFR and NFE vs those that do not")
gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(avg_cds_len = mean(cds_length),
            median_cds_len = median(cds_length),
            avg_gene_len = mean(gene_length),
            median_gene_len = median(gene_length),
            avg_n_exons = mean(num_coding_exons),
            median_n_exons = median(num_coding_exons),
            avg_n_sites = mean(n_sites),
            median_n_sites = median(n_sites),
            avg_maf_af = mean(max_af),
            median_maf_af = median(max_af)
            )

# Next step to add to script:
# Create further subcategories of differences 







gnomad_lof %>% 
  ggplot(aes(x = gene_length, y = num_coding_exons)) +
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = cds_length, y = num_coding_exons)) +
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = gene_length, y= n_sites)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = cds_length, y = n_sites)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = num_coding_exons, y= oe_lof)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = num_coding_exons, y= oe_lof_upper)) + 
  geom_point(alpha = 0.1) + 
  theme_bw()

gnomad_lof %>% 
  ggplot(aes(x = p_nfe, y = p_afr, col = diff)) + 
#  geom_abline(slope =1, intercept = 0, color = "firebrick")+  
 # geom_line(data = above_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
#  geom_line(data = below_df_line, aes(x = p_nfe, y = predicted), color = "blue", size = 0.5, linetype = "dashed") +
  geom_point() + 
  theme_bw() +
  labs(title = "AFR vs NFE (p)") 



