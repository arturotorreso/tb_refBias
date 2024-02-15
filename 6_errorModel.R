rm(list = ls())

library(vcfR)
library(seqinr)
library(foreach)
library(parallel)
library(doSNOW)

dir_prefix = '/'

meta_file = paste(dir_prefix, 'lineages_simplified.txt', sep = '')
sample_file = paste(dir_prefix, 'assemblies_final.txt', sep = '')
ref_file = paste(dir_prefix, 'Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa', sep = '')
genes_file = paste(dir_prefix, 'h37rv_genes2.bed', sep = '')


metadata = read.table(meta_file, header = T, stringsAsFactors = F)
refseq = read.fasta(ref_file)$Chromosome
refgenes = read.delim(genes_file, header = F, stringsAsFactors = F, col.names = c('st', 'end', 'gene', 'gene_id'))


sample_df = read.table(sample_file, stringsAsFactors = F, col.names = 'sample')
metadata = metadata[metadata$sample %in% sample_df$sample,]
metadata = metadata[grep('1|2|3|4', metadata$lineage),]


all_lineages = read.delim("assembly_lineages.txt", header = F, stringsAsFactors = F)[,1:2]
colnames(all_lineages) = c('id', 'lineage')
all_lineages$lineage[all_lineages$lineage == ''] = 'M.microti'

all_lineages = all_lineages[all_lineages$id %in% metadata$sample[metadata$lineage == 'lineage4'],]
all_lineages$lineage = sapply(all_lineages$lineage, function(x) paste(strsplit(strsplit(x, ';')[[1]][length(strsplit(x, ';')[[1]])],'[.]')[[1]][1:2], collapse = '.'))


  
## 1. Load sample and null distribution

# Load sample distribution
load(paste(dir_prefix, '/1_mapping/mapping_compVCF.RData', sep = ''))
# load(paste(dir_prefix, '/2_assembly_mapping/ass_mapp_compVCF.RData', sep = ''))

# Load the null distribution
mapping_null = readRDS('~/Documents/PhD/TB/ref_genome_analysis/assembly/1_mapping/null_dist.mapping.rds')
# assembly_null = readRDS('~/Desktop/null_dist.assembly.rds')

gene_df_all = gene_df_all[gene_df_all$lineage == 'lineage4',]
# mapping_null = mapping_null[mapping_null$lineage == 'lineage4',]

rownames(all_lineages) = all_lineages$id

gene_df_all$lineage = all_lineages[gene_df_all$sample,'lineage']


## 2. Error model (error vs null dis)

out_all$sample[is.na(out_all$sample)] = out_all$ref[is.na(out_all$sample)]
out_all$assembly[is.na(out_all$assembly)] = out_all$ref[is.na(out_all$assembly)]


# Sample distribution
gene_model = gene_df_all[,c('gene', 'st', 'end', 'length', 'sample', 'lineage', 'rep')]
gene_model$error = (gene_df_all$fn + gene_df_all$fp) / gene_df_all$length



n.cores = 6
cl <- makeCluster(n.cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = length(unique(gene_model$gene)), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


gene_glm = foreach (gene=unique(gene_model$gene), .combine='rbind', .options.snow = opts) %:%
           # foreach (lineage=unique(metadata$lineage), .combine='rbind') %dopar% {
             foreach (lineage=unique(gene_model$lineage), .combine='rbind') %dopar% {
    gene_model2 = gene_model[gene_model$gene == gene & gene_model$lineage == lineage,]
    null_dist = mapping_null[mapping_null$lineage == strsplit(lineage, '[.]')[[1]][1],]
    # null_dist = mapping_null
    colnames(null_dist)[1] = 'gene'
    null_dist$gene = paste(gene, '_null',sep = '')
    gene_model2 = rbind(gene_model2[,c('gene', 'lineage', 'error')], null_dist[,c('gene', 'lineage', 'error')])
    gene_model2$gene = factor(gene_model2$gene, levels = c(paste(gene, '_null', sep = ''), gene))
    
    m1 = glm(error ~ gene, data = gene_model2)
    
    # m1_ci = suppressMessages(confint(m1))
    
    m1_null = as.numeric(m1$coefficients[1])
    m1_sample = m1$coefficients[1] + m1$coefficients[2]
    
    # ci_null = m1_ci[1,]
    # ci_sample = m1$coefficients[1] + m1_ci[2,]
    
    p = summary(m1)$coefficients[2,4]
    return(c(gene, lineage, m1_null, m1_sample, p))
    # return(c(gene, lineage, m1_null, ci_null, m1_sample, ci_sample, p))
}


close(pb)
stopCluster(cl)

gene_glm = as.data.frame(gene_glm)
colnames(gene_glm) = c('gene', 'lineage', 'null', 'null_low', 'null_upp', 'sample', 'sample_low', 'sample_upp', 'p_value')


# write.table(gene_glm, 'geneModel_mapping.txt',
#             row.names = F, col.names = F, quote = F, sep = '\t')
write.table(gene_glm, 'geneModel_assembly_all.txt',
            row.names = F, col.names = F, quote = F, sep = '\t')

