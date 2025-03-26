

library(enrichR)


############### Get gene sets to test enrichment for ####################

gnomad_p_diff = data.table::fread("gnomad_lof_metrics_delta_p_0.1_na_rm.csv")

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
                   dplyr::slice_sample(n = length(input_genes)) %>% 
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
         "GO_Molecular_Function_2025",
         "GO_Biological_Process_2025",
         "Rare_Diseases_GeneRIF_Gene_Lists", #associated with rare disease on pubmed search - not text-mining
         "GWAS_Catalog_2023", 
         "OMIM_Disease",
         "MGI_Mammalian_Phenotype_Level_4_2024",
         "ClinVar",
         "GTEx_Tissues_V8_2023"
         )

# writeLines(input_genes, "gene_symbols.txt")
# writeLines(background_genes, "background_gene_symbols.txt")

####################### Perform enrichment ########################


get_nfe_diff_genes("afr", c("amr", "fin"))

enrichment = enrichr(
  genes = input_genes,
  databases = dbs,
  background = background_genes,
  include_overlap = TRUE,
  sleepTime = 10
)


enrichment$Chromosome_Location %>% 
  ggplot(aes(x = Adjusted.P.value)) + 
  geom_histogram() + 
  xlim(0,1) +
  theme_bw()

enrichment$`Epigenomics_Roadmap_HM_ChIP-seq` %>% head
enrichment$GTEx_Tissues_V8_2023
enrichment$GO_Molecular_Function_2025 %>%  dplyr::arrange(Adjusted.P.value) %>% head()

enrichment$Human_Phenotype_Ontology  %>% dplyr::arrange(Adjusted.P.value) %>% head()

enrichment$GTEx_Tissues_V8_2023 %>% dplyr::arrange(Adjusted.P.value) %>% head()

enrichment$GTEx_Tissue_Expression_Up %>% dplyr::arrange(Adjusted.P.value) %>% head()

enrichment$OMIM_Disease %>% dplyr::arrange(Adjusted.P.value) %>% head()
enrichment$Chromosome_Location %>% dplyr::arrange(Adjusted.P.value) %>% head()

enrichment$GWAS_Catalog_2023 %>% dplyr::arrange(Adjusted.P.value) %>%  head()

enrichment$UK_Biobank_GWAS_v1 %>% dplyr::arrange(Adjusted.P.value) %>% head()
