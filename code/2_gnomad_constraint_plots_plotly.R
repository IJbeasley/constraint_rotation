
# Following up from the previous script, we will now use the plotly package to make the scatter plots of
# p vs p_nfe (across populations) interactive.
#  This is used to help use explore the features of genes where there are obvious population differences in p.

# Load the data
# using gnomad v2.1.1 
gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

library(ggplot2)
library(dplyr)

############ getting overall distribution of important values to compare to the subset of genes with differences in p ###############

{

print("\n All populations p distribution")
summary(gnomad_lof$p)

print("\n Transcript level distribution")
print("Transcript level from Gencode")
summary(as.factor(gnomad_lof$transcript_level))

print("\n CDs length distribution")
summary(gnomad_lof$cds_length)

print("\n Gene length distribution")
summary(gnomad_lof$gene_length)

print("\n N sites distribution")
summary(gnomad_lof$n_sites)

print("\n Brain expression distribution")
print("Brain expression from GTEx")
summary(as.factor(gnomad_lof$brain_expression))

}

############# Scatter plot of p vs p_nfe #############

lof_plot_afr_nfe = gnomad_lof %>% 
  ggplot(aes(x=p_afr, 
             y=p_nfe,
              # adding all these extra aes arguments below
              #  so we can se these gene features when we hover over the points
              # in the interactive plotly plot
             name = gene, 
             cds_length= cds_length,
             gene_length = gene_length,
             n_distinct_lof_sites = n_sites,
             num_coding_exons = num_coding_exons,
             gene_type = gene_type,
             transcript_level = transcript_level,
             transcript_type = transcript_type,
             brain_expression = brain_expression,
             chrom = chromosome
  )) + 
  geom_abline(slope = 1, intercept = 0, col = "red") + 
  geom_point() + 
  theme_bw() +
  labs(title = "AFR vs NFE (p)") 

plotly::ggplotly(lof_plot_afr_nfe)

###### scatter plot of p_fin vs p_nfe ######

lof_plot_afr_nfe= gnomad_lof %>% 
  ggplot(aes(x=p_fin, 
             y=p_nfe,
             name = gene, 
             cds_length= cds_length,
             gene_length = gene_length,
             n_distinct_lof_sites = n_sites,
             num_coding_exons = num_coding_exons,
             gene_type = gene_type,
             transcript_level = transcript_level,
             transcript_type = transcript_type,
             brain_expression = brain_expression,
             chrom = chromosome
  )) + 
  geom_abline(slope = 1, intercept = 0, col = "red") + 
  geom_point() + 
  theme_bw() +
  labs(title = "FIN vs  NFE (p)")

plotly::ggplotly(lof_plot_afr_nfe)

############ scatter plot of p_amr vs p_nfe ############

lof_plot_afr_nfe = gnomad_lof %>% 
  ggplot(aes(x=p_amr, 
             y=p_nfe,
             name = gene, 
             cds_length= cds_length,
             gene_length = gene_length,
             n_distinct_lof_sites = n_sites,
             num_coding_exons = num_coding_exons,
             gene_type = gene_type,
             transcript_level = transcript_level,
             transcript_type = transcript_type,
             brain_expression = brain_expression,
             chrom = chromosome
  )) + 
  geom_abline(slope = 1, intercept = 0, col = "red") + 
  geom_point() + 
  theme_bw() +
  labs(title = "AMR vs NFE (p)")

plotly::ggplotly(lof_plot_afr_nfe)
