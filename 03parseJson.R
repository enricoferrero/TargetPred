library(jsonlite)

dat <- stream_in(gzfile("../data/16.04_association_data.json.gz"))
dat <- cbind.data.frame(dat$target$id, dat$target$gene_info$symbol, dat$disease$id, dat$disease$efo_info$label, dat$is_direct, dat$association_score$overall, dat$association_score$datatypes)
names(dat) <- c("EnsemblId", "Symbol", "OntologyId", "Label", "Is.direct", "overall", "literature", "rna_expression", "genetic_association", "somatic_mutation", "known_drug", "animal_model", "affected_pathway")
write.csv(dat, gzfile("../data/open_targets_association_data.csv.gz"), row.names=FALSE)
