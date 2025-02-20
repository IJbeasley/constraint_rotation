

gnomad_lof = data.table::fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")

library(ggplot2)
library(dplyr)

lof_plot_afr_nfe = gnomad_lof %>% 
  ggplot(aes(x=p_afr, 
             y=p_nfe, # adding all these aes for the plotly graph interactive
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


summary(as.factor(gnomad_lof$transcript_level))
summary(gnomad_lof$cds_length)
summary(gnomad_lof$n_sites)

plotly::ggplotly(lof_plot_afr_nfe)
