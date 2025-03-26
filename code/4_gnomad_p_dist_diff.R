
# And
# Quantify comparisons in gene features for 
# genes that have different p between AFR and NFE groups
# vs genes that do not


{
############################ Load the data ########################
# using gnomad v2.1.1

gnomad_lof = data.table::fread("output/p_diff/gnomad_lof_metrics_delta_p_0.1_na_rm.csv")

library(ggplot2)
library(dplyr)

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







##################### Overlap between constraint differences ############

{
gnomad_lof = gnomad_lof %>% 
  group_by(afr_nfe_p_diff)

message("Median of gene length features, 
        for genes with different constraint 
        between AFR and NFE vs those that do not")

print(
  gnomad_lof %>% 
    summarise(median_cds_len = median(cds_length),
              median_gene_len = median(gene_length),
              median_n_exons = median(num_coding_exons)
    )
)
}

gnomad_lof = gnomad_lof %>% 
  group_by(afr_nfe_p_diff,
           amr_nfe_p_diff,
           fin_nfe_p_diff)

gnomad_lof %>% 
  summarise(n = n())



############################### Compare genes with different p values between populations ###############################

# Do they have different gene features to those that do not differ between AFR and NFE?


{
message("\n Number of genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  # group_by(diff) %>% 
  summarise(n = n())
)

message("\n Distribution of p_afr and p_nfe for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  # group_by(diff) %>%
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
  # group_by(diff) %>%
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
#  group_by(diff) %>% 
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
  # group_by(diff) %>% 
  summarise(n = n(), 
                    percent_testis_gene = mean(is_testis_gene)
                   )
)

# summarise_gnomad_lof = gnomad_lof %>% 
#   group_by(afr_nfe_p_diff,
#               amr_nfe_p_diff,
#               fin_nfe_p_diff,
#              is_testis_gene
#            ) %>% 
#     summarise(n = n())
# 
# 
# 
# summarise_gnomad_lof = summarise_gnomad_lof %>% 
#   dplyr::filter(!is.na(afr_nfe_p_diff)) %>% 
#   mutate(is_testis_gene = ifelse(is_testis_gene == 1, "testis_gene", "not_testis_gene"),
#          afr_nfe_p_diff = ifelse(afr_nfe_p_diff == 1, "afr_nfe_diff", "afr_nfe_same"),
#          amr_nfe_p_diff = ifelse(amr_nfe_p_diff == 1, "amr_nfe_diff", "amr_nfe_same"),
#          fin_nfe_p_diff = ifelse(fin_nfe_p_diff == 1, "fin_nfe_diff", "fin_nfe_same")
#          ) %>% 
#   tidyr::pivot_wider(values_from = "n",
#                      names_from = c(afr_nfe_p_diff,
#                                      amr_nfe_p_diff,
#                                      fin_nfe_p_diff
#                                      )
#   )
# 
# # Convert to matrix
# contingency_table <- as.matrix(summarise_gnomad_lof[, -1])


diff = gnomad_lof %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(p)) %>% 
  dplyr::group_by(afr_nfe_p_diff, is_testis_gene) %>% 
  summarise(n = n())
  

contingency_table <- matrix(
  diff$n,  # Fill values row-wise
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    afr_nfe_p_diff = c("0", "1"),  # Row labels
    is_testis_gene = c("0", "1")   # Column labels
  )
)

fisher.test(contingency_table)

# just below 0.05 significance ... 
print(
  pbinom(q = floor(0.425 * 73), 
       size = 73, 
       prob = 0.328, 
       lower.tail = F
       ) 
)
      
print("\n Distribution of max MAF for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
  # group_by(diff) %>% 
  summarise(avg_maf_af = mean(max_af),
                    median_maf_af = median(max_af)
  )
)

message("\n Distribution of gene length features, for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof %>% 
#  group_by(diff) %>% 
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
  # group_by(diff) %>% 
  summarise(
            avg_n_sites = mean(n_sites),
            median_n_sites = median(n_sites)
            )
)
  
} 

########################## Check patterns for just highly constrained genes ############# 
# Next step to add to script:
# Create further subcategories of differences 

# Define big threshold as close to median p 
big = 0.2

gnomad_lof = gnomad_lof %>% dplyr::ungroup()

gnomad_lof_filt = gnomad_lof %>% 
                  dplyr::filter(p_afr < big)

gnomad_lof_filt = gnomad_lof_filt %>% 
  group_by(afr_nfe_p_diff)
           #amr_nfe_p_diff,
           #fin_nfe_p_diff)

# gnomad_lof  = gnomad_lof  %>% 
#                         dplyr::mutate(
#                                               unconstrained_p_afr =ifelse(p_afr >= big, 1, 0),
#                                               unconstrained_p_nfe = ifelse(p_nfe >= big, 1, 0),
#                                               
#                                               ) %>% 
#                       dplyr:: group_by(diff, unconstrained_p_afr, unconstrained_p_nfe)                        

gnomad_lof_filt %>% 
  summarise(n = n())

message("\n Distribution of gene length features, for genes that differ between AFR and NFE vs those that do not")
print(
  gnomad_lof_filt %>% 
    #  group_by(diff) %>% 
    summarise(median_gene_len = median(gene_length),
              median_n_exons = median(num_coding_exons),
              median_n_sites = median(n_sites, na.rm = T)
    )
)




