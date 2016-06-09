library(jsonlite)

dat <- stream_in(gzfile("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v3.0/16.04_association_data.json.gz"))
dat <- cbind.data.frame(dat$target$id, dat$disease$id, dat$is_direct, dat$association_score$datatypes)
names(dat) <- c("EnsemblId", "OntologyId", "Is.direct", "literature", "rna_expression", "genetic_association", "somatic_mutation", "known_drug", "animal_model", "affected_pathway")
write.csv(dat, gzfile("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v3.0/public.matrix.csv.gz"), quote=FALSE, row.names=FALSE)
