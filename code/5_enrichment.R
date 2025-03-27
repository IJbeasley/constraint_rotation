

library(enrichR)
# library(speedr)
# speedr::set_server("https://maayanlab.cloud/enrichrapi")
# speedr::list_libraries()
# 
# 
# url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GWAS_Catalog_2023"
# 
# download.file(url, "data/enrichr/GWAS_Catalog_2023.gmt")
# 
# res <-  speedr::read_gmt("data/enrichr/GWAS_Catalog_2023.gmt")
# 
# result <- speedr::enrich("GWAS_Catalog_2023", input_genes)

############### Get gene sets to test enrichment for ####################

gnomad_p_diff = data.table::fread("output/p_diff/gnomad_lof_metrics_delta_p_0.1_na_rm.csv")

# get genes that have different estimated constraint to nfe for given pop
# and same / similar estimated constraint to nfe for other pops
get_nfe_diff_genes = function(diff_pop, #population that should be different to nfe
                              same_pop # population that should be same to nfe
                              ) {


  
# get different genes
input_genes <- gnomad_p_diff %>% 
               dplyr::filter(!!sym(paste0(diff_pop,"_nfe_p_diff")) == 1) %>% 
               dplyr::filter(dplyr::if_all(paste0(same_pop,"_nfe_p_diff"),
                                           ~.x == 0
                                           )
                             ) %>% 
              dplyr::pull(gene)

# get background genes - one per value of p for diff pop
background_genes <- gnomad_p_diff %>% 
                    dplyr::filter(!!sym(paste0(diff_pop,"_nfe_p_diff")) == 0) %>% 
                    dplyr::group_by(!!sym(paste0("p_", diff_pop))) %>% 
                    dplyr::slice_sample(n = 1) 


background_genes = background_genes %>%
                   dplyr::ungroup() %>% 
                   #dplyr::slice_sample(n = length(input_genes)) %>% 
                   dplyr::pull(gene)



assign("input_genes", input_genes, envir = .GlobalEnv)
assign("background_genes", background_genes, envir = .GlobalEnv)

message("\n Number of genes to test enrichment")
print(length(input_genes))
message("\n Number of background genes")
print(length(background_genes))

} 


############## What gene set enrichment to perform ###################

dbs <- c(
        #"Chromosome_Location",
        # "Human_Phenotype_Ontology",
         # "GO_Molecular_Function_2025",
         # "GO_Biological_Process_2025",
         "Rare_Diseases_GeneRIF_Gene_Lists", #associated with rare disease on pubmed search - not text-mining
         "GWAS_Catalog_2023", 
         "OMIM_Disease",
         "MGI_Mammalian_Phenotype_Level_4_2024",
       #  "ClinVar",
         "GTEx_Tissues_V8_2023"
         )

################ Make + save gene sets ################
{


get_nfe_diff_genes("afr", c("amr", "fin"))

writeLines(input_genes, 
           "output/enrichr/afr_diff_nfe_gene_symbols.txt"
           )

writeLines(background_genes, 
           "output/enrichr/afr_diff_nfe_background_gene_symbols.txt"
           )



get_nfe_diff_genes("amr", c("afr", "fin"))

writeLines(input_genes, 
           "output/enrichr/amr_diff_nfe_gene_symbols.txt"
)

writeLines(background_genes, 
           "output/enrichr/amr_diff_nfe_background_gene_symbols.txt"
)



get_nfe_diff_genes("fin", c("amr", "afr"))

writeLines(input_genes, 
           "output/enrichr/fin_diff_nfe_gene_symbols.txt"
)

writeLines(background_genes, 
           "output/enrichr/fin_diff_nfe_background_gene_symbols.txt"
)


}

####################### Perform enrichment ########################

pops = c("afr", "amr", "fin")
pops = "amr"

pop_colors <- c(
  "AFR" = "#941494", 
  "AMR" = "#ED1E24", 
  "FIN" = "#002060"
)

for(pop in pops){
  
  
input_genes = readLines(paste0("output/enrichr/",
                               pop,
                               "_diff_nfe_gene_symbols.txt")
                        )


enrichment = enrichr(
  genes = input_genes,
  databases = dbs,
#  background = background_genes,
  include_overlap = TRUE,
  sleepTime = 10
)





for(gene_set_test in dbs){
  
  if(grepl("GTEx", gene_set_test)){
    
    n_to_plot = 15
  } else {
    
    n_to_plot = 10
  }
  
  sorted_data =  enrichment[[gene_set_test]] %>% 
                 arrange(Adjusted.P.value, P.value) %>% 
                 dplyr::slice_head(n = n_to_plot)
  
  gene_set_test_name = gsub("_", " ", gene_set_test)
  gene_set_test_name = paste(toupper(pop), " - ", gene_set_test_name)

plot = ggplot(sorted_data, 
       aes(x = factor(Term, levels = rev(sorted_data$Term)),
             #Term, #reorder(Term, P.value, decreasing = T), 
           y = Adjusted.P.value)) +
  geom_col(fill = "white", colour = pop_colors[toupper(pop)]) +
  coord_flip() +  # Flip axes for better readability
  geom_text(aes(label = Term, 
                y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
                ),
            position = position_stack(vjust = 0.5)
            ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),  # Remove major y-axis grid lines
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
        ) +
  labs(y = "Adjusted P-value", 
       x = "",
       title = gene_set_test_name)


print(plot)

ggsave(paste0("presentation_figs/", pop, "_", gene_set_test, ".png"), 
       width = 17, 
       height = 14, 
       units = "cm"
) 
  
}

}


sorted_data = enrichment$GTEx_Tissues_V8_2023 %>%
  dplyr::select(Term, P.value, Adjusted.P.value, Genes) %>%
  #dplyr::arrange(Genes) %>%
  dplyr::group_by(Genes, Adjusted.P.value) %>%
  #dplyr::arrange(P.value, Adjusted.P.value) %>%
  dplyr::slice_sample(n = 1)

sorted_data = sorted_data %>% 
              dplyr::ungroup() %>% 
              dplyr::arrange(P.value, Adjusted.P.value) %>% 
              dplyr::slice_head(n = 5)

enrichment$GTEx_Tissues_V8_2023 %>% 
              dplyr::filter(Genes == "P2RX5;TTC24;MADCAM1")

tissue_lab = c("Small Intestine", "Spleen", "Kidney - Cortex", "Vagina", "Adipose")

sorted_data$Tissue = tissue_lab

sorted_data$n_demo_groups = c(1,8,1,1,8)

sorted_data = sorted_data %>% 
              dplyr::mutate(label = paste0(Tissue, "\n (", Genes, ")"))

enrichment$GTEx_Tissues_V8_2023 %>% 
  dplyr::filter(Genes == "P2RX5;TTC24;NOX5")

# Male 20-29, 30-39, 40-49, 60 - 69 Female 20-29, 40-49, 50-59, 60 - 69
# grouped_data = enrichment$GTEx_Tissues_V8_2023 %>% 
#   dplyr::select(Term, P.value, Adjusted.P.value, Genes) %>% 
#   #dplyr::arrange(Genes) %>%  
#   dplyr::group_by(Genes, Adjusted.P.value, P.value) %>% 
#   #dplyr::arrange(P.value, Adjusted.P.value) %>% 
#   dplyr::slice_sample(n = 1)

plot = ggplot(sorted_data, 
              aes(x = factor(Genes, levels = rev(sorted_data$Genes)),
                  #Term, #reorder(Term, P.value, decreasing = T), 
                  y = Adjusted.P.value)) +
  geom_col(fill = "white", colour = pop_colors[toupper(pop)]) +
  coord_flip() +  # Flip axes for better readability
  geom_text(aes(label = label,
                y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  ),
  position = position_stack(vjust = 0.5)
  ) +
  # geom_text(aes(label = Genes,
  #               y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  # ),
  # position = position_stack(vjust = 1)
  # ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),  # Remove major y-axis grid lines
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
  ) +
  labs(y = "Adjusted P-value", 
       x = "",
       title = gene_set_test_name)

print(plot)

pop = "afr"
ggsave(paste0("presentation_figs/", pop, "_", gene_set_test, ".png"), 
       width = 17, 
       height = 14, 
       units = "cm"
) 




enrichment$GWAS_Catalog_2023 %>% arrange(Adjusted.P.value, P.value) %>% View()

################## AMR


pop = "amr"
sorted_data = enrichment$GTEx_Tissues_V8_2023 %>%
  dplyr::select(Term, P.value, Adjusted.P.value, Genes) %>%
  #dplyr::arrange(Genes) %>%
  dplyr::group_by(Genes, Adjusted.P.value) %>%
  #dplyr::arrange(P.value, Adjusted.P.value) %>%
  dplyr::slice_sample(n = 1)

sorted_data = sorted_data %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(Adjusted.P.value, P.value) %>% 
  dplyr::slice_head(n = 5)


enrichment$GTEx_Tissues_V8_2023 %>% 
  dplyr::filter(Genes == "EFCAB13;PKHD1L1")

enrichment$GTEx_Tissues_V8_2023 %>% 
  dplyr::filter(Genes == "COL6A5")

enrichment$GTEx_Tissues_V8_2023 %>% 
  dplyr::filter(Genes == "SLC5A9;ENPP7")

sorted_data$n_demo_groups = c(1,1,10,6,1)

tissue_lab = c("Thyroid","Small Intestine", "Lung", "Adrenal Gland", "Liver")

sorted_data$Tissue = tissue_lab

sorted_data = sorted_data %>% 
  dplyr::mutate(label = paste0(Tissue, "\n (", Genes, ")"))



plot = ggplot(sorted_data, 
              aes(x = factor(Genes, levels = rev(sorted_data$Genes)),
                  #Term, #reorder(Term, P.value, decreasing = T), 
                  y = Adjusted.P.value)) +
  geom_col(fill = "white", colour = pop_colors[toupper(pop)]) +
  coord_flip() +  # Flip axes for better readability
  geom_text(aes(label = label,
                y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  ),
  position = position_stack(vjust = 0.5)
  ) +
  # geom_text(aes(label = Genes,
  #               y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  # ),
  # position = position_stack(vjust = 1)
  # ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),  # Remove major y-axis grid lines
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
  ) +
  labs(y = "Adjusted P-value", 
       x = "",
       title = gene_set_test_name)

print(plot)

pop = "amr"
ggsave(paste0("presentation_figs/", pop, "_", gene_set_test, ".png"), 
       width = 17, 
       height = 14, 
       units = "cm"
) 


############################### FIN



pop = "fin"
sorted_data = enrichment$GTEx_Tissues_V8_2023 %>%
  dplyr::select(Term, P.value, Adjusted.P.value, Genes) %>%
  #dplyr::arrange(Genes) %>%
  dplyr::group_by(Genes, Adjusted.P.value) %>%
  #dplyr::arrange(P.value, Adjusted.P.value) %>%
  dplyr::slice_sample(n = 1)

sorted_data = sorted_data %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(Adjusted.P.value, P.value) %>% 
  dplyr::slice_head(n = 5)

sorted_data$n_demo_groups = c(3,4,10)

tissue_lab = c("Minor Salivary Gland", "Spleen", "Pituitary")

sorted_data$Tissue = tissue_lab

sorted_data = sorted_data %>% 
  dplyr::mutate(label = paste0(Tissue, "\n (", Genes, ")"))


plot = ggplot(sorted_data, 
              aes(x = factor(Genes, levels = rev(sorted_data$Genes)),
                  #Term, #reorder(Term, P.value, decreasing = T), 
                  y = Adjusted.P.value)) +
  geom_col(fill = "white", colour = pop_colors[toupper(pop)]) +
  coord_flip() +  # Flip axes for better readability
  geom_text(aes(label = label,
                y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  ),
  position = position_stack(vjust = 0.5)
  ) +
  # geom_text(aes(label = Genes,
  #               y = Adjusted.P.value #+ 0.02 * max(Adjusted.P.value)
  # ),
  # position = position_stack(vjust = 1)
  # ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),  # Remove major y-axis grid lines
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)
  ) +
  labs(y = "Adjusted P-value", 
       x = "",
       title = gene_set_test_name)

print(plot)

pop = "fin"
ggsave(paste0("presentation_figs/", pop, "_", gene_set_test, ".png"), 
       width = 17, 
       height = 14, 
       units = "cm"
) 

