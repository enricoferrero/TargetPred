options(stringsAsFactors=FALSE)
library(biomaRt)
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")

## gene_target
gtres <- read.delim("../data/docstore.gene_target.results.txt")
gtres$symbol <- sub("^.+GENE\\$_?([[:alnum:]]+).+$", "\\1", gtres$sentence)

# annotate genes
genes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="external_gene_name", values=gtres$symbol, mart=ensembl)
gtres <- merge(gtres, genes, by.x="symbol", by.y="external_gene_name", all.x=TRUE, all.y=FALSE)

# export
gtres <- gtres[c("ensembl_gene_id", "symbol", "docId")]
names(gtres) <- c("ensembl", "symbol", "pubmed")
gtres <- unique(na.omit(gtres))
write.csv(gtres, "../data/docstore.gene_target.annotated.csv", quote=FALSE, row.names=FALSE)
