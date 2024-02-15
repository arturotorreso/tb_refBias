rm(list = ls())

library(vcfR)
library(seqinr)
library(doSNOW)
library(foreach)

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


## 1. Load Raw mapping bias data 
load(paste(dir_prefix, '/1_mapping/mapping_compVCF.RData', sep = ''))
mapping_gene_df_all = gene_df_all
mapping_out_all = out_all

mapping_out_all$sample[is.na(mapping_out_all$sample)] = mapping_out_all$ref[is.na(mapping_out_all$sample)]
mapping_out_all$assembly[is.na(mapping_out_all$assembly)] = mapping_out_all$ref[is.na(mapping_out_all$assembly)]


load(paste(dir_prefix, '/2_assembly_mapping/ass_mapp_compVCF.RData', sep = ''))
assembly_gene_df_all = gene_df_all
assembly_out_all = out_all

assembly_out_all$sample[is.na(assembly_out_all$sample)] = assembly_out_all$ref[is.na(assembly_out_all$sample)]
assembly_out_all$assembly[is.na(assembly_out_all$assembly)] = assembly_out_all$ref[is.na(assembly_out_all$assembly)]


refgenes$length = refgenes$end - refgenes$st


out_all = assembly_out_all
# out_all = mapping_out_all

################################################################################



mapping_null = c()
for (lineage in unique(metadata$lineage)){
  print(lineage)
  samples = metadata$sample[metadata$lineage == lineage]
  samples = if (length(samples)>=40) sample(samples, 40) else samples
  
  # Null distribution
  len = range(refgenes$length)
  len = seq(round(len[1],-1), round(len[2],-1), 1000)
  
  n.cores = 6
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "FORK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  cat(paste('Parellel environment: ', foreach::getDoParRegistered(), '\n',sep = ''))
  cat(paste('Number of nodes: ', foreach::getDoParWorkers(), '\n', sep = ''))
  
  null_lineage = foreach (l=len, .combine='rbind') %:%
  foreach (i = 1:100, .combine = 'rbind') %dopar% {
      print(l)
      st = sample(seq(1, 4411532 - l), 1)
      end = st + l
      null_dist_sample = c()
      for (s in samples){
        for (n in sample(1:10, 3)) {
          out_all2 = out_all[out_all$id == s & out_all$rep == n,]
          
          snps = out_all2[out_all2$pos >= st & out_all2$pos <= end,]
          
          fn = ifelse(nrow(snps) > 0, nrow(snps[is.na(snps$sample),]), 0)
          fp = ifelse(nrow(snps) > 0, nrow(snps[is.na(snps$assembly),]), 0)
          
          error = (fp + fn)/l
          
          null_dist_sample = rbind(null_dist_sample, c(paste(l, '_null', sep = ''), lineage, fn, fp, error))
        }
      }
      null_dist_sample
    }

  parallel::stopCluster(cl = my.cluster)
  
  
  mapping_null = rbind(mapping_null, null_lineage)
  
}


mapping_null = as.data.frame(mapping_null)
colnames(mapping_null) = c('length', 'lineage', 'fn', 'fp', 'error')
mapping_null$fn = as.numeric(mapping_null$fn)
mapping_null$fp = as.numeric(mapping_null$fp)
mapping_null$error = as.numeric(mapping_null$error)


# saveRDS(mapping_null, '~/Desktop/null_dist.mapping.rds')
saveRDS(mapping_null, '~/Desktop/null_dist.assembly.rds')

