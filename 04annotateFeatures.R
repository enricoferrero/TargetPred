library(biomaRt)
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")

fv <- readRDS(file.path("../data/fv.rds"))
fv <- fv$data
fv$name <- sub("feat.", "", fv$name)
fv$name <- sub("GO_", "GO:", fv$name)
fv$name <- sub("R_HSA_", "R-HSA-", fv$name)

interpro <- fv[grep("^IPR", fv$name), "name"]
interpro <- getBM(attributes=c("interpro", "interpro_description"), filters="interpro", values=interpro, mart=ensembl)
names(interpro) <- c("name", "description")
interpro <- na.omit(merge(fv, interpro, all=TRUE))

go_id <- fv[grep("^GO", fv$name), "name"]
go_id <- getBM(attributes=c("go_id", "name_1006"), filters="go_id", values=go_id, mart=ensembl)
names(go_id) <- c("name", "description")
go_id <- na.omit(merge(fv, go_id, all=TRUE))

reactann <- read.delim("http://www.reactome.org/download/current/ReactomePathways.txt", header=FALSE)
reactome <- fv[grep("^R-HSA", fv$name), "name"]
reactome <- unique(merge(reactome, reactann[1:2], by=1, all.x=TRUE, all.y=FALSE))
names(reactome) <- c("name", "description")
reactome <- na.omit(merge(fv, reactome, all=TRUE))

efoann <- read.csv("/GWD/bioinfo/projects/bix-analysis-stv/data/CTTV/v2.0/matrix.csv.gz")
efoann <- efoann[c("OntologyId", "Label")]
efo <- fv[grep("^genetic_association|^somatic_mutation|^rna_expression|^affected_pathway|^animal_model", fv$name), "name"]
efo <- as.data.frame(do.call(rbind, strsplit(efo, "\\.")))
efo <- unique(merge(efo, efoann, by.x=2, by.y=1, all.x=TRUE, all.y=FALSE))
efo$name <- paste0(efo[[2]], ".", efo[[1]])
efo$description <- paste0(efo[[3]], " ", "(", efo[[2]], ")")
efo <- na.omit(merge(fv, efo[c("name", "description")], all=TRUE))

others <- fv[grep("^IPR|^GO|R-HSA|^genetic_association|^somatic_mutation|^rna_expression|^affected_pathway|^animal_model", fv$name, invert=TRUE), ]
others$description <- others$name

feats <- rbind(interpro, go_id, reactome, efo, others)
feats <- feats[order(feats[3], decreasing=TRUE), ]
write.csv(feats, file.path("../data/SelectedFeatures.csv"), quote=TRUE, row.names=FALSE)
