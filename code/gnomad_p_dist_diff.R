

# Quantify comparisons in gene features for 
# genes that have different p between AFR and NFE groups
# vs genes that do not

# Next step to add to script:
# Create further subcategories of differences 



gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(diff = ifelse(p_afr >= above_formula | p_afr <= below_formula, 1, 0)) %>% 
  dplyr::mutate(diff = as.factor(diff))

gnomad_lof %>% 
  group_by(diff) %>% 
  summarise(n = n())

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

gtex_med_tpm = data.table::fread("data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct")

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
  labs(x= "PCA Components")

library(plotly)

plot = pca$rotation %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Tissue") %>%  
  dplyr::mutate(Testis = ifelse(Tissue == "Testis", 1, 0)) %>% 
  dplyr::mutate(Testis = as.factor(Testis)) %>% 
  ggplot(aes(x=PC1, y=PC6, col = Testis, Tissue = Tissue)) + 
  geom_point() + 
  theme_bw() 

ggplotly(plot)

testis_expr <- gtex_med_tpm %>%
               dplyr::pull(Testis) #%>%

pc_scores <- as.data.frame(pca$x)

# Compute correlations with Testis expression
pc_cor <- cor(pc_scores, testis_expr, use = "pairwise.complete.obs", method = "spearman")

pc_cor_df <- data.frame(PC = colnames(pc_scores), Correlation = pc_cor)

# Sort by absolute correlation
pc_cor_df <- pc_cor_df %>% arrange(desc(abs(Correlation)))

# Print top PCs
print(pc_cor_df[1:5,])

testis_genes = pca$x %>% 
  as.data.frame() %>% 
  dplyr::select(PC6) %>% 
  dplyr::filter(PC6 < -0.5) %>% 
  row.names()
  #head()

  
  
  
pca$x %>% 
  as.data.frame()  %>% 
  row.names() %>% 
  unique() %>% 
  length()


pca$x %>% 
  as.data.frame()  %>% 
  row.names() %>% 
  length()



# Remove the dot and number after it
testis_genes <- sub("\\.\\d+$", "", testis_genes)


gnomad_lof = gnomad_lof %>% 
  dplyr::mutate(is_testis_gene = ifelse(gene_id %in% testis_genes, 1, 0)) #%>% 
  # group_by(is_testis_gene) %>% 
  # summarise(n = n())

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

