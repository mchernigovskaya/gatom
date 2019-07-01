refine_datatable <- function(dt) {
    as.data.table(as.data.frame(dt),
                  keep.rownames = !is.numeric(attr(dt, "row.names")))
}


prepareData <- function(gene.de.raw, top=12000, k.positive.scores=200) {
    gene.de.unique <- unique(gene.de.raw[order(pval)], by="symbol")
    gene.de.unique <- gene.de.unique[order(-baseMean)][1:top] # TODO head
    gene.de.unique[, score := -log(pval)]
    threshold <- gene.de.unique[order(-score)][k.positive.scores, score]
    gene.de.unique[, score := score - threshold]
}


prepareMet <- function(met.de.raw, k.positive.scores=200) {
    HMDB2metabolite <- fread("./inst/reactome/hmdb2kegg.tsv")
    setkey(HMDB2metabolite, HMDB)
    setkey(met.de.raw, ID)
    met.de.raw <- HMDB2metabolite[met.de.raw, nomatch=0]
    met.de.unique <- unique(met.de.raw[order(pval)], by="KEGG")
    met.de.unique[, score := -log(pval)]
    threshold <- met.de.unique[order(-score)][k.positive.scores, score]
    met.de.unique[, score := score - threshold]
}



coveredReactions <- function(gene.de.unique, reactome) {
    catalyst <- reactome[grepl("catalystActivity", types),]
    reaction_coverage <- catalyst[,.(isCovered=any(symbol %in% gene.de.unique[, symbol])), by=reaction_id] ##! unique
    reaction_coverage[isCovered == TRUE, reaction_id]
}


prepareReactome <- function(reactome, gene.de.unique, mask=TRUE) {
    uniprot <- reactome[db_type=="UniProt"]
    # UNIPROT -> ENTREZ -> SYMBOL
    uniprot[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                              AnnotationDbi::mapIds(org.Mm.eg.db, uniprot[,cross_ref], 'ENTREZID', 'UNIPROT'),
                                              'SYMBOL', 'ENTREZID')]
    compound <- reactome[db_type=="COMPOUND"]
    compound[, symbol := compound[,cross_ref]] # whatever
    reactome <- rbindlist(list(uniprot, compound))
    reactome <- reactome[sapply(symbol, function(x)!is.null(x))]
    reactome[, symbol := unlist(symbol)]
    setkeyv(reactome, c("reaction_id", "symbol"))
    reactome <- unique(reactome)
    reactome <- reactome[reaction_id %in% coveredReactions(gene.de.unique, reactome), ]
    reactome[, types_short:=substr(types, 1, 1)]
    if (mask) {
        reactome <- maskMets(reactome)
        reactome <- maskGenes(reactome)
    }
    reactome
}
