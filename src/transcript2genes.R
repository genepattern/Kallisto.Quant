print('R-->Running the Rscript now')

destfile='/build/for_sleuth/transcript2genes_dataframe.RDA'
if(!file.exists(destfile)){
print('R-->About to setup biomaRt')
library("biomaRt")
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ext_gene = external_gene_name)
print('R-->biomaRt has been setup')
save(t2g,file="transcript2genes_dataframe.RDA")
print("t2g has been saved to 'transcript2genes_dataframe.RDA' ")
print(head(t2g))
} else {
print('R-->About to load a previous biomaRt run')
load(destfile)
print(head(t2g))
}

print('R-->About to call sleuth')
library(ggplot2) # req'd by sleuth
library(dplyr, warn.conflicts = FALSE) # req'd by sleuth; masks some functions, but this doesn't affect us
library("sleuth")
s2c = do.call(rbind.data.frame, list(c('temp_name_1', 'temp_name_2', 'RNASeq_quant/')))
colnames(s2c) <- c('sample', 'condition', 'path')
s2c[] <- lapply(s2c, as.character)
so <- sleuth_prep(s2c, ~condition, target_mapping=t2g, aggregation_column='ext_gene',gene_mode=TRUE)
print('R-->Done with sleuth')

write.csv(so$obs_norm[, c(2,4)], 'RNASeq_quant/gene_expression.csv', quote=FALSE, row.names=TRUE)
