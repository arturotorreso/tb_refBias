rm(list = ls())

library(vcfR)
library(seqinr)
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


######## 1. Load Raw mapping bias data ########

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


###########################################################
###########################################################



######## 2. Mapping bias per lineage ########

mapping_model = read.delim('/Users/arturo/Documents/PhD/TB/ref_genome_analysis/assembly/1_mapping/geneModel_mapping.txt', header = F, stringsAsFactors = F)
colnames(mapping_model) = c('gene', 'lineage','null', 'null_low', 'null_upp', 'sample', 'sample_low', 'sample_upp', 'p_value')

assembly_model = read.delim('/Users/arturo/Documents/PhD/TB/ref_genome_analysis/assembly/2_assembly_mapping/geneModel_assembly.txt', header = F, stringsAsFactors = F)
colnames(assembly_model) = c('gene', 'lineage','null', 'null_low', 'null_upp', 'sample', 'sample_low', 'sample_upp', 'p_value')


gene_ann = read.delim('/Users/arturo/Documents/PhD/TB/ref_genome_analysis/metadata/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt', header = T,
                      stringsAsFactors = F)



mapping_model = mapping_model[mapping_model$lineage != 'lineage4',]
assembly_model = assembly_model[assembly_model$lineage != 'lineage4',]

# mapping_model$lineage = factor(mapping_model$lineage, levels = c('lineage4', 'lineage1', 'lineage2', 'lineage3'))
mapping_model$lineage = factor(mapping_model$lineage, levels = c('lineage4.9', unique(mapping_model$lineage)[!unique(mapping_model$lineage) %in% 'lineage4.9']))
assembly_model$lineage = factor(assembly_model$lineage, levels = c('lineage4.9', unique(assembly_model$lineage)[!unique(assembly_model$lineage) %in% 'lineage4.9']))


# summary(m1 <- glm(sample~lineage + gene, data = mapping_model[mapping_model$gene %in% sample(unique(mapping_model$gene),500),]))
summary(m1 <- glm(sample~lineage, data = mapping_model))

lin_error = cbind(names(coef(m1)), coef(m1), confint(m1))
lin_error = as.data.frame(lin_error)
colnames(lin_error) = c('lineage', 'fit', 'low', 'upp')
rownames(lin_error) = NULL
lin_error$lineage[lin_error$lineage == '(Intercept)'] = 'lineagelineage4.9'
lin_error$lineage = sub('lineage', '', lin_error$lineage)
class(lin_error$fit) = class(lin_error$low) = class(lin_error$upp)= 'numeric'
lin_error$fit[-1] = lin_error$fit[-1] + lin_error$fit[1]
lin_error$low[-1] = lin_error$low[-1] + lin_error$fit[1]
lin_error$upp[-1] = lin_error$upp[-1] + lin_error$fit[1]
lin_error = lin_error[order(lin_error$lineage),]


lin_error$low[lin_error$low < 0] = 0


summary(m2 <- glm(sample~lineage, data = assembly_model))

lin_error2 = cbind(names(coef(m2)), coef(m2), confint(m2))
lin_error2 = as.data.frame(lin_error2)
colnames(lin_error2) = c('lineage', 'fit', 'low', 'upp')
rownames(lin_error2) = NULL
lin_error2$lineage[lin_error2$lineage == '(Intercept)'] = 'lineagelineage4.9'
lin_error2$lineage = sub('lineage', '', lin_error2$lineage)
class(lin_error2$fit) = class(lin_error2$low) = class(lin_error2$upp)= 'numeric'
lin_error2$fit[-1] = lin_error2$fit[-1] + lin_error2$fit[1]
lin_error2$low[-1] = lin_error2$low[-1] + lin_error2$fit[1]
lin_error2$upp[-1] = lin_error2$upp[-1] + lin_error2$fit[1]
lin_error2 = lin_error2[order(lin_error2$lineage),]


lin_error2$low[lin_error2$low < 0] = 0



pdf('lin_error.pdf', width = 8, height = 3.5)

par(mar = c(4,7,3,3))
plot(lin_error$fit, type = 'n', axes = F, ann=F, ylim = c(min(lin_error$low), max(lin_error$upp)))
box(lwd = 2)
axis(1, at = 1:nrow(lin_error), labels = toupper(sub('ineage', '', lin_error$lineage)))
axis(2, las = 2)
mtext('Error rate', 2, 5, font = 2, cex = 1.3)

points(1:11 - 0.10, lin_error$fit, pch = 19, col = 'grey', cex = 1.3)
arrows(1:11 - 0.10, lin_error$low, 1:11 - 0.10, col = 'grey', lin_error$upp, code = 3, angle = 90, length = 0.05, lwd= 1.3)

points(1:11 + 0.10, lin_error2$fit, pch = 19, col = 'black', cex = 1.3)
arrows(1:11 + 0.10, lin_error2$low, 1:11 + 0.10, col = 'black', lin_error2$upp, code = 3, angle = 90, length = 0.05, lwd= 1.3)

text(x = 2.2, y = 0.000185, 'Mapping', col = 'grey', font = 2, adj = 0, cex = 0.8)
text(x = 2.2, y = 0.000168, 'Assembly + Mapping', col = 'black', font = 2, adj = 0, cex = 0.8)


dev.off()


write.table(mapping_model, 'mapping_table.tsv', sep = '\t', row.names = F, quote = F)





######## 3. Mapping bias vs nsnps ########

nsamples = 100
samples = unique(mapping_out_all$id)
samples = sample(samples,nsamples)

genome_size = 4411532


n.cores = 5
cl <- makeCluster(n.cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsamples, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


out = foreach (i = 1:100, .combine = 'rbind', .options.snow = opts) %dopar% {
  s = samples[i]

  n = sample(1:10,1)
  
  wsize = 5
  st = sample(1:1000, 1)
  end = st + wsize
  entry = mapping_out_all[mapping_out_all$id == s & mapping_out_all$rep == n,]
  lin = unique(mapping_out_all$lineage[mapping_out_all$id == s])
  wout = c()
  while (end < (genome_size + wsize)){
    wentry = entry[entry$pos >= st & entry$pos <= end,]
    wentry = wentry[wentry$assembly != wentry$ref | wentry$sample != wentry$ref,]
    
    if (nrow(wentry) == 0){
      st = end + 1
      end = st + wsize
      next()
    }
    
    nsnps = nrow(wentry)
    fn = nrow(wentry[wentry$sample == wentry$ref & wentry$assembly != wentry$ref,])
    fp = nrow(wentry[wentry$sample != wentry$ref & wentry$assembly == wentry$ref,])
    
    
    error = fp + fn
    wout = rbind(wout, c(s, n, st, end, nsnps, fn, fp, error))
    
    st = end + 1
    end = st + wsize
  }
  wout
}

close(pb)
stopCluster(cl) 


out = as.data.frame(out)
colnames(out) = c('id', 'rep', 'st', 'end', 'nsnps', 'fn', 'fp', 'error')

out$error = as.numeric(out$error)
out$nsnps = as.numeric(out$nsnps)


out$perror = out$error/out$nsnps
plot(out$nsnps, out$perror)

df = out[out$id %in% unique(out$id),c('id','nsnps', 'perror')]
colnames(df) = c('id','x', 'y')

mlog <- nls(y ~ SSlogis(x, phi1, phi2, phi3),
              data = df)

library(nlme)
mmlog = nlme::nlme(y ~ SSlogis(x, phi1, phi2, phi3),
            random= phi1+phi2+phi3~1|id,
            fixed = phi1+phi2+phi3~1,
            data = df,
            start = c(coef(mlog)[1], coef(mlog)[2], coef(mlog)[3]),
            control = lmeControl(opt = "optim", maxIter=1000))



interval = sapply(unique(df$id), function(id) SSlogis(seq(0,10,0.01), coef(mmlog)[id,1], coef(mmlog)[id,2], coef(mmlog)[id,3]))


pdf('error_snp_window.pdf', height = 5, width = 5, useDingbats = F)
par(mar = c(6,6,3,3))
plot(seq(0,10,0.01), SSlogis(seq(0,10,0.01), mmlog$coefficients$fixed[1], mmlog$coefficients$fixed[2], mmlog$coefficients$fixed[3]),
     type = 'n', lwd = 2, xlab= '', ylab = '', axes=FALSE,ann=FALSE, ylim = c(0,1))
box(lwd = 2)
axis(1, cex.axis = 1.2)
axis(2, las = 2, cex.axis = 1.2)
mtext('Number of SNPs in genomic window (10bp)', 1, 3, font = 2, cex = 1.1)
mtext('Proportion of Errors', 2, 3.5, font = 2, cex = 1.2)

polygon(x = c(seq(0,10,0.01), rev(seq(0,10,0.01))), y = c(apply(interval, 1, function(x) quantile(x, 0.025)), rev(apply(interval, 1, function(x) quantile(x, 0.975)))), col = adjustcolor('steelblue', 0.5),border = F)
lines(seq(0,10,0.01), SSlogis(seq(0,10,0.01), mmlog$coefficients$fixed[1], mmlog$coefficients$fixed[2], mmlog$coefficients$fixed[3]), lwd = 2)
dev.off()



## Gene family


gene_ann = read.delim('metadata/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt', header = T,
                      stringsAsFactors = F)

mapping_model$gene[mapping_model$gene == 'Rv2418c'] = 'octT'
mapping_gene_df_all$gene[mapping_gene_df_all$gene == 'Rv2418c'] = 'octT'

assembly_model$gene[assembly_model$gene == 'Rv2418c'] = 'octT'
assembly_gene_df_all$gene[assembly_gene_df_all$gene == 'Rv2418c'] = 'octT'


mapping_model = merge(x = mapping_model, y = gene_ann[,c('Name', 'Functional_Category')], by.x = 'gene', by.y = 'Name', all.x = T, all.y = F)


sapply(split(mapping_model$sample, mapping_model$Functional_Category), mean) 

mapping_model$Functional_Category = factor(mapping_model$Functional_Category,
                                           levels = c('stable RNAs',
                                                      unique(mapping_model$Functional_Category)[!unique(mapping_model$Functional_Category) %in% 'stable RNAs']))

summary(m1 <- glm(sample ~ Functional_Category, data = mapping_model))






mapping_gene_df_all$error = (mapping_gene_df_all$fn + mapping_gene_df_all$fp) / mapping_gene_df_all$length
assembly_gene_df_all$error = (assembly_gene_df_all$fn + assembly_gene_df_all$fp) / assembly_gene_df_all$length

mapping_gene_df_all = merge(x = mapping_gene_df_all, y = gene_ann[,c('Name', 'Functional_Category')], by.x = 'gene', by.y = 'Name', all.x = T, all.y = F)
assembly_gene_df_all = merge(x = assembly_gene_df_all, y = gene_ann[,c('Name', 'Functional_Category')], by.x = 'gene', by.y = 'Name', all.x = T, all.y = F)


mapping_gene_df_all$lineage = factor(mapping_gene_df_all$lineage, c('lineage4', 'lineage1','lineage2','lineage3'))
mapping_gene_df_all$Functional_Category = factor(mapping_gene_df_all$Functional_Category,
                                           levels = c('stable RNAs',
                                                      unique(mapping_gene_df_all$Functional_Category)[!unique(mapping_gene_df_all$Functional_Category) %in% 'stable RNAs']))


summary(m3 <- lm(error ~ lineage*Functional_Category, data = mapping_gene_df_all))


pred_model = data.frame('lineage' = rep(unique(mapping_gene_df_all$lineage), length(unique(mapping_gene_df_all$Functional_Category))), 
           'Functional_Category' = rep(unique(mapping_gene_df_all$Functional_Category), length(unique(mapping_gene_df_all$lineage))))
pred_model$lineage = as.character(pred_model$lineage)
pred_model = pred_model[order(pred_model$Functional_Category, pred_model$lineage),]
pred_model = cbind(pred_model, predict(m3,newdata = pred_model, interval="confidence", level=0.95))






assembly_gene_df_all$lineage = factor(assembly_gene_df_all$lineage, c('lineage4', 'lineage1','lineage2','lineage3'))
assembly_gene_df_all$Functional_Category = factor(assembly_gene_df_all$Functional_Category,
                                                 levels = c('stable RNAs',
                                                            unique(assembly_gene_df_all$Functional_Category)[!unique(assembly_gene_df_all$Functional_Category) %in% 'stable RNAs']))

summary(m4 <- lm(error ~ lineage*Functional_Category, data = assembly_gene_df_all))

pred_model2 = data.frame('lineage' = rep(unique(assembly_gene_df_all$lineage), length(unique(assembly_gene_df_all$Functional_Category))), 
                        'Functional_Category' = rep(unique(assembly_gene_df_all$Functional_Category), length(unique(assembly_gene_df_all$lineage))))
pred_model2$lineage = as.character(pred_model2$lineage)
pred_model2 = pred_model2[order(pred_model2$Functional_Category, pred_model2$lineage),]
pred_model2 = cbind(pred_model2, predict(m4,newdata = pred_model2, interval="confidence", level=0.95))


pred_model$lwr[pred_model$lwr < 0] = 0
pred_model2$lwr[pred_model2$lwr < 0] = 0





pred_model_log = pred_model
pred_model2_log = pred_model2

pred_model_log$fit = log10(pred_model_log$fit)
pred_model2_log$fit = log10(pred_model2_log$fit)

pred_model_log$lwr = log10(pred_model_log$lwr)
pred_model2_log$lwr = log10(pred_model2_log$lwr)

pred_model_log$upr = log10(pred_model_log$upr)
pred_model2_log$upr = log10(pred_model2_log$upr)




###################
###  ALL PLOTS  ###
###################


pdf('~/Documents/PhD/TB/ref_genome_analysis/manuscript/figures/fig1_new2.pdf', height = 7, width = 9.5, useDingbats = F)

layout(matrix(c(1,4,4,4,4,
                4,4,4,4,4,
                2,2,2,3,9,
                5,5,5,9,9,
                6,6,6,9,9,
                7,7,7,9,9,
                8,8,8,9,9
),nrow = 7, byrow = T),
heights = c(.3,2,.3,1,1,1,1))

par(oma = c(6,2,2,2))
# layout.show(9)

for (l in c('a','b','c')){
  if (l == 'a'){par(mar = c(0,0,0,0))} else if (l == 'b') {par(mar = c(0,0,1,0))} else {par(mar = c(0,4,1,0))}
  plot.new()
  text(0,1,l, font = 2, cex = 2.5, xpd = NA)
}


## Fig A
par(mar = c(4,8,2,3)) 
plot(lin_error$fit, type = 'n', axes = F, ann=F, ylim = c(min(lin_error$low), max(lin_error$upp)))
box(lwd = 1.4)
axis(1, at = 1:nrow(lin_error), labels = toupper(sub('ineage', '', lin_error$lineage)), cex.axis = 1.3)
axis(2, las = 2, cex.axis = 1.3)
mtext('Error rate', 2, 6, font = 2, cex = 1.2, xpd = NA)

points(1:nrow(lin_error) - 0.10, lin_error$fit, pch = 19, col = 'grey', cex = 1.3)
arrows(1:nrow(lin_error) - 0.10, lin_error$low, 1:nrow(lin_error) - 0.10, col = 'grey', lin_error$upp, code = 3, angle = 90, length = 0.05, lwd= 1.3)

points(1:nrow(lin_error) + 0.10, lin_error2$fit, pch = 19, col = 'black', cex = 1.3)
arrows(1:nrow(lin_error) + 0.10, lin_error2$low, 1:nrow(lin_error) + 0.10, col = 'black', lin_error2$upp, code = 3, angle = 90, length = 0.05, lwd= 1.3)


text(x = 2.2, y = 0.000193, 'Mapping', col = 'grey', font = 2, adj = 0, cex = 1)
text(x = 2.2, y = 0.000170, 'Assembly + Mapping', col = 'black', font = 2, adj = 0, cex = 1)



## Fig B
par(mar = c(1,7,0,0))

x = unique(pred_model_log$Functional_Category)
y = rep(10, length(x))

ylabs = c('0.0001', '0.0005', '0.001')
ypos = log10(c(0.0001, 0.0005, 0.001))

xlabs = c("Stable RNAs", "Conserved hypotheticals", 
  "Virulence, detoxification, adaptation", "Metabolism and respiration", 
  "Lipid metabolism", "Cell wall and cell processes", "Information pathways", 
  "Regulatory proteins", "PE/PPE", "Insertion seqs and phages", 
  "Unknown")

for (i in 1:4){
  lineage = unique(pred_model_log$lineage)[i]
  plot(x,y,type = 'n', axes = F, ann=F,
       ylim = c(min(c(pred_model_log$lwr,pred_model2_log$lwr)),
                max(c(pred_model_log$upr, pred_model2_log$upr)))) 
  box(lwd = 2)
  axis(2, las = 2, at = ypos, labels = ylabs)
  
  text(paste('Lineage ', i, sep = ''), x=1, y=-3, adj = 0.5, font = 2, cex = 1)
  
  if (i == 1){
    legend(x = 3, y = -3, legend = c('Mapping', 'Assembly + Mapping'), pch = 19, col = c('steelblue1', 'steelblue4'), bty = "n", yjust = 0.5)
  }
  
  
  if (i == 4){
    axis(1, at = 1:length(unique(pred_model2_log$Functional_Category)),
         labels = rep('', length(unique(pred_model2_log$Functional_Category))))
    text(x = 1:length(unique(pred_model2_log$Functional_Category)),
         y = par("usr")[3] - 0.23,
         labels = xlabs,
         xpd = NA,
         srt = 20,
         adj = 0.965,
         cex = 1, font = 2)
  }
  
  for (i in 1:11){
    lines(c(i, i),
          c(pred_model_log$fit[pred_model_log$lineage == lineage][i],
            pred_model2_log$fit[pred_model2_log$lineage == lineage][i]),
          col = 'steelblue')
  }
  
  points(1:11, pred_model_log$fit[pred_model_log$lineage == lineage], pch = 19, col = 'steelblue1', cex = 1.3)
  points(1:11, pred_model2_log$fit[pred_model2_log$lineage == lineage], pch = 19, col = 'steelblue4', cex = 1.3)
  
  
}

mtext('Error rate', 2, -2, font = 2, cex = 1, outer = TRUE, adj = 0.25)


## Fig C

par(mar = c(13,10,4,3))
plot(seq(0,10,0.01), SSlogis(seq(0,10,0.01), mmlog$coefficients$fixed[1], mmlog$coefficients$fixed[2], mmlog$coefficients$fixed[3]),
     type = 'n', lwd = 2, xlab= '', ylab = '', axes=FALSE,ann=FALSE, ylim = c(0,1))
box(lwd = 1.4)
grid(lwd = 1.3) 
axis(1, cex.axis = 1.2)
axis(2, las = 2, cex.axis = 1.2)
mtext('Number of SNPs\nin genomic window (10bp)', 1, 4.5, font = 2, cex = 0.8)
mtext('Proportion of Errors', 2, 3.5, font = 2, cex = 0.8)

polygon(x = c(seq(0,10,0.01), rev(seq(0,10,0.01))), y = c(apply(interval, 1, function(x) quantile(x, 0.025)), rev(apply(interval, 1, function(x) quantile(x, 0.975)))), col = adjustcolor('steelblue', 0.5),border = F)
lines(seq(0,10,0.01), SSlogis(seq(0,10,0.01), mmlog$coefficients$fixed[1], mmlog$coefficients$fixed[2], mmlog$coefficients$fixed[3]), lwd = 2)

dev.off()

