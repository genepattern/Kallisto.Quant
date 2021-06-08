#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

print('R-->Running the Rscript now')
print('R-->Arguments read are')
print(args)

library(rhdf5)
## Assumption: Kallisto output file is RNASeq_quant/
## Original post was deleted: https://www.biostars.org/p/419605/
oldids<-h5read('RNASeq_quant/abundance.h5','/aux/ids')
print(head(oldids))
print('-----------------')
newids <- gsub("\\|.*", "", oldids)
print(head(newids))
h5write(newids, 'RNASeq_quant/abundance.h5', "/aux/ids")

print('R-->About to load a previous transcript2genes file')
t2g <- read.csv(file =args[1],colClasses='character')
print(head(t2g))

print('R-->About to call sleuth')

library(ggplot2) # req'd by sleuth
library(dplyr, warn.conflicts = FALSE) # req'd by sleuth; masks some functions, but this doesn't affect us
library("sleuth")
s2c = do.call(rbind.data.frame, list(c('gp_sample', 'gp_sample', 'RNASeq_quant/')))
colnames(s2c) <- c('sample', 'condition', 'path')
s2c[] <- lapply(s2c, as.character)
so <- sleuth_prep(s2c, ~condition, target_mapping=t2g, aggregation_column='gene_name',gene_mode = TRUE)

print("Privew of normalized data:")
print(head(so$obs_norm))

print("Privew of transcript level (raw) data:")
print(head(so$obs_raw))

print("--------------------")
print("Privew of gene level (raw) estimated counts:")
counts = sleuth_to_matrix(so, 'obs_norm', 'est_counts')
print(head(counts))

# print("--------------------")
# tpm = sleuth_to_matrix(so, 'obs_norm', 'tpm')
# print("sleuth_to_matrix(so, 'obs_norm', 'tpm')")
# print(head(tpm))

write.csv(so$obs_norm[, c(2,4)], paste0('RNASeq_quant/',args[2],'_normalized_gene_level.csv'), quote=FALSE, row.names=FALSE)
# write.csv(so$obs_norm, paste0('RNASeq_quant/',args[2],'_full_gene_level.csv'), quote=FALSE, row.names=FALSE)
write.csv(so$obs_raw, paste0('RNASeq_quant/',args[2],'_raw_transcript_level.csv'), quote=FALSE, row.names=FALSE)
write.csv(counts, paste0('RNASeq_quant/',args[2],'_raw_estimated_counts.csv'), quote=FALSE, row.names=TRUE)
