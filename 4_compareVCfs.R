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


### Compare VCFs

n.cores = 6
cl <- makeCluster(n.cores)
registerDoSNOW(cl)
pb <- utils::txtProgressBar(max = length(metadata$sample), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


out_all = foreach (s=metadata$sample, .combine='rbind', .options.snow = opts) %dopar% {
  vcf_assembly_file = paste(dir_prefix, 'assemblies_h37rv/vcf/',s,'.refH37rv.vcf.gz', sep = '')
  vcf_assembly = vcfR::read.vcfR(vcf_assembly_file, verbose = F)
  vcf_assembly = vcf_assembly[!grepl('INDEL', vcf_assembly@fix[,'INFO'])]
  
  lineage = metadata[metadata$sample == s,'lineage']
  out_sample = c()
  for (n in 1:10) {
    vcf_sample_file = paste(dir_prefix, '2_assembly_mapping/vcf_onlyContigs/',s,'.rep',n,'.sfilt.vcf.gz', sep = '')
    vcf_sample = vcfR::read.vcfR(vcf_sample_file, verbose = F)
    
    vcf_sample = vcf_sample[vcf_sample@fix[,'FILTER'] == 'PASS']
    vcf_sample = vcf_sample[!grepl('INDEL', vcf_sample@fix[,'INFO'])]
    
    pos = as.numeric(unique(c(vcf_assembly@fix[,'POS'], vcf_sample@fix[,'POS'])))
    if (length(pos) == 0){next()}
    out = c()
    for (p in pos){
      ref = toupper(refseq[p])
      
      assembly_entry = vcf_assembly@fix[vcf_assembly@fix[,"POS"] == p]
      sample_entry = vcf_sample@fix[vcf_sample@fix[,"POS"] == p]
      
      if (length(assembly_entry) > 0){
        assembly_dp = as.numeric(strsplit(grep('DP=', strsplit(assembly_entry[8], ';')[[1]], value = T), '=')[[1]][2])
        if (assembly_dp > 1) {next()}
      }
      
      assembly_ref = ref
      assembly_alt = ref
      if (length(assembly_entry) > 0){
        assembly_ad = as.numeric(strsplit(strsplit(vcf_assembly@gt[vcf_assembly@fix[,"POS"] == p][2],':')[[1]][which(strsplit(vcf_assembly@gt[vcf_assembly@fix[,"POS"] == p][1],':')[[1]]=='AD')], ',')[[1]])
        assembly_gt = unique(as.numeric(strsplit(strsplit(vcf_assembly@gt[vcf_assembly@fix[,"POS"] == p][2],':')[[1]][which(strsplit(vcf_assembly@gt[vcf_assembly@fix[,"POS"] == p][1],':')[[1]]=='GT')], '/')[[1]])) + 1
        if (length(assembly_gt) > 1){
          assembly_gt = which(assembly_ad == max(assembly_ad))
        }
        assembly_ref = assembly_entry[4]
        assembly_alt = paste(c(assembly_ref, strsplit(assembly_entry[5], ',')[[1]])[assembly_gt], collapse = '/')
      }
      
      sample_ref = ref
      sample_alt = ref
      if (length(sample_entry) > 0){
        sample_ad = as.numeric(strsplit(strsplit(vcf_sample@gt[vcf_sample@fix[,"POS"] == p][2],':')[[1]][which(strsplit(vcf_sample@gt[vcf_sample@fix[,"POS"] == p][1],':')[[1]]=='AD')], ',')[[1]])
        sample_gt = unique(as.numeric(strsplit(strsplit(vcf_sample@gt[vcf_sample@fix[,"POS"] == p][2],':')[[1]][which(strsplit(vcf_sample@gt[vcf_sample@fix[,"POS"] == p][1],':')[[1]]=='GT')], '/')[[1]])) + 1
        if (length(sample_gt) > 1){
          sample_gt = which(sample_ad == max(sample_ad))
        }
        sample_ref = sample_entry[4]
        sample_alt = paste(c(sample_ref, strsplit(sample_entry[5], ',')[[1]])[sample_gt], collapse = '/')
      }
      
      out = rbind(out, c(s, n, lineage, p, ref, assembly_alt, sample_alt, isTRUE(assembly_alt==sample_alt)))
    }
    
    out_sample = rbind(out_sample, out)
  }
  
  return(out_sample)
}


out_all = as.data.frame(out_all)
colnames(out_all) = c('id', 'rep', 'lineage', 'pos', 'ref', 'assembly', 'sample', 'shared')
out_all$pos = as.integer(out_all$pos)





n.cores = 6
cl <- makeCluster(n.cores)
registerDoSNOW(cl)
pb <- utils::txtProgressBar(max = length(metadata$sample), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


gene_df_all = foreach (s=metadata$sample, .combine='rbind', .options.snow = opts) %dopar% {
  lineage = metadata[metadata$sample == s,'lineage']
  gene_df_n = c()
  for (n in 1:10) {
    out = out_all[out_all$id == s & out_all$rep == n,]
    
  
    gene_df = refgenes[,c('gene', 'st', 'end')]
    gene_df$length = gene_df$end - gene_df$st 
    gene_df$sample = s
    gene_df$lineage = lineage
    gene_df$rep = n
    gene_df$fn = 0
    gene_df$fp = 0
    # gene_df$tn = NA
    # gene_df$tp = NA
    
    if (nrow(out) == 0){
      next()
    }
    
    for (i in 1:nrow(refgenes)){
      gene = gene_df$gene[i]
      st = gene_df$st[i]
      end = gene_df$end[i]
      len = gene_df$length[i]
      
      snps = out[out$pos >= st & out$pos <= end,]
      
      gene_df$fn[i] = ifelse(nrow(snps) > 0, nrow(snps[snps$sample == snps$ref & snps$assembly != snps$ref,]), 0)
      gene_df$fp[i] = ifelse(nrow(snps) > 0, nrow(snps[snps$sample != snps$ref & snps$assembly == snps$ref,]), 0)
      
    }
    
    gene_df_n = rbind(gene_df_n, gene_df)
  }
  
  return(gene_df_n)
}


# save(out_all, gene_df_all, file = paste(dir_prefix, 'ref_genome_analysis/assembly/1_mapping/mapping_compVCF.RData', sep = ''))
save(out_all, gene_df_all, file = paste(dir_prefix, 'ref_genome_analysis/assembly/2_assembly_mapping/ass_mapp_compVCF.RData', sep = ''))



