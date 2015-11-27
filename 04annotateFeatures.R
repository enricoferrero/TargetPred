library(biomaRt)
#ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="www.ensembl.org")

# separate small molecules and antibodies
for (agenttype in c("small_molecule", "antibody")) {

    fv <- readRDS(file.path(paste0("../data/fv.", agenttype, ".rds")))
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
    reactome <- unique(merge(reactome, reactann[1:2], by=1))
    names(reactome) <- c("name", "description")
    reactome <- na.omit(merge(fv, reactome, all=TRUE))

    others <- fv[grep("^(?!IPR)", fv$name, perl=TRUE), ]
    others <- others[grep("^(?!GO)", others$name, perl=TRUE), ]
    others <- others[grep("^(?!R-HSA)", others$name, perl=TRUE), ]
    others$description <- others$name

    feats <- rbind(interpro, go_id, reactome, others)
    feats <- feats[order(feats[3], decreasing=TRUE), ]
    write.csv(feats, file.path(paste0("../data/SelectedFeatures.", agenttype, ".csv")), quote=TRUE, row.names=FALSE)

}









