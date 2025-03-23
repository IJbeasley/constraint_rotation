

library(enrichR)

gnomad_p_diff = data.table::fread("gnomad_lof_metrics_delta_p_0.1_na_rm.csv")

input_genes <- gnomad_p_diff %>% 
               dplyr::filter(afr_nfe_p_diff == 1 & amr_nfe_p_diff == 0 & fin_nfe_p_diff == 0) %>% 
               dplyr::pull(gene)

gnomad_p_diff %>% 
  dplyr::filter(afr_nfe_p_diff == 1) %>% 
  dplyr::pull(gene_id)

background_genes <- gnomad_p_diff %>% 
                    dplyr::filter(afr_nfe_p_diff == 0) %>% 
                    #dplyr::slice_sample(n = 100) %>% 
                    dplyr::pull(gene)


message("\n Number of genes to test enrichment")
print(length(input_genes))
message("\n Number of background genes")
print(length(background_genes))

dbs <- c("Chromosome_Location",
         "Human_Phenotype_Ontology",
         "GO_Molecular_Function_2025",
         "GO_Biological_Process_2025",
         "Rare_Diseases_GeneRIF_Gene_Lists", #associated with rare disease on pubmed search - not text-mining
         "GWAS_Catalog_2023", 
         "OMIM_Disease",
         "MGI_Mammalian_Phenotype_Level_4_2024",
         "Epigenomics_Roadmap_HM_ChIP-seq",
         "TF-LOF Expression from GEO",
         "ClinVar",
         "GTEx_Tissues_V8_2023",
         "Human_Gene_Atlas"
         )

writeLines(input_genes, "gene_symbols.txt")
writeLines(background_genes, "background_gene_symbols.txt")

dbs = "ChEA_2022"

input_genes = "MAML3"

enrichment = enrichr(
  genes = input_genes,
  databases = dbs,
  #background = background_genes,
  include_overlap = TRUE,
  sleepTime = 1
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
