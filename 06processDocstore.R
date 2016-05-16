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
write.table(gtres, "../data/docstore.gene_target.annotated.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## gene_target_indication
# process results
gtires <- read.delim("../data/docstore.gene_target_indication.results.txt")
gtires$symbol <- sub("^.+GENE\\$_?([[:alnum:]]+).+$", "\\1", gtires$sentence)
gtires$mesh <- sub("^.+INDICATION\\$_?([[:alnum:]]+).+$", "\\1", gtires$sentence)

# annotate genes
genes <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters="external_gene_name", values=gtires$symbol, mart=ensembl)
gtires <- merge(gtires, genes, by.x="symbol", by.y="external_gene_name", all.x=TRUE, all.y=FALSE)

# annotate diseases
diseases <- read.delim("../data/mesh2efo.mappings.txt")
diseases$MESH <- sub("^.+/", "", diseases$MESH)
diseases$EFO <- sub("^.+/", "", diseases$EFO)
gtires <- merge(gtires, diseases, by.x="mesh", by.y="MESH", all.x=TRUE, all.y=FALSE)

# export
gtires <- gtires[c("ensembl_gene_id", "symbol", "EFO", "mesh", "docId")]
names(gtires) <- c("ensembl", "symbol", "efo", "mesh", "pubmed")
gtires <- unique(na.omit(gtires))
write.table(gtires, "../data/docstore.gene_target_indication.annotated.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
