
solution2module <- function(solution, gene.de.unique, met.de.unique) {
    m <- graph.data.frame(solution, directed=FALSE)
    print(m)
    mets <- metabolitesInfo(V(m)$name)
    mets_names <- mets[match(V(m)$name, mets[, metabolite]), metabolite_name]
    V(m)$label <- ifelse(!is.na(mets_names), mets_names, V(m)$name)

    #V(m)$logPval <- gene.de.unique[match(V(m)$name, gene.de.unique[,symbol]), log(pval)]

    V(m)$logPval <- gene.de.unique[match(V(m)$name, gene.de.unique[,symbol]), log(pval_gene)]
    V(m)$logPvalPhosph <- gene.de.unique[match(V(m)$name, gene.de.unique[,symbol]), log(pval_phosph)]

    V(m)$log2FC <- gene.de.unique[match(V(m)$name, gene.de.unique[,symbol]), log2FC]
    V(m)$url <- mets[match(V(m)$name, mets[, metabolite]), metabolite_url]

    if (!is.null(met.de.unique)) {
        mets_pval <- met.de.unique[match(V(m)$name, met.de.unique[,KEGG]), log(pval)]
        V(m)$logPval <- ifelse(!is.na(mets_pval), mets_pval, V(m)$logPval)
        mets_log2FC <- met.de.unique[match(V(m)$name, met.de.unique[,KEGG]), log2FC]
        V(m)$log2FC <- ifelse(!is.na(mets_log2FC), mets_log2FC, V(m)$log2FC)
    }
    m
}


modulePostProcessing_raw <- function(m, reactome, interesting_nodes=NA) {
    # add reactions url
    V(m)$url <- ifelse(startsWith(V(m)$name, "R-"), sprintf("https://reactome.org/content/detail/%s", V(m)$name), V(m)$url)

    # add labels to edges
    elabel <- unique(reactome[,.(reaction_id, symbol, types_short)])
    elabel <- elabel[, label:=list(typ=paste0(types_short, collapse=";")), by=.(reaction_id, symbol)]
    elabel[, edge:=paste(reaction_id, symbol, sep=" ")]
    elabel[, rev_edge:=paste(symbol, reaction_id, sep=" ")]
    elist <- get.edgelist(m)
    elist <- paste(elist[,1], elist[,2])
    #E(m)$name <- elabel[match(elist, elabel[, edge]), label]
    E(m)$label <- elabel[match(elist, elabel[, edge]), label]
    rev_names <- elabel[match(elist, elabel[, rev_edge]), label] # WTF? Why (gene, reactione) if I add (r, v)
    E(m)$label <- ifelse(!is.na(rev_names), rev_names, E(m)$label)

    # add reaction names
    reaction_names <- data.table(read_csv("inst/reactome/mus_reaction_names.csv"))
    reaction_names[, unique_name:=paste(reaction_name, reaction_id)]
    reaction_names <- reaction_names[match(V(m)$name, reaction_names[, reaction_id]), unique_name]
    V(m)$name <- ifelse(!is.na(reaction_names), reaction_names, V(m)$name)

    # add shape
    is_met <- match(V(m)$name, reactome[, cross_ref])
    V(m)$nodeType <- ifelse(!is.na(is_met), "circle", "square")

    # reactions is circle
    V(m)$nodeType <- ifelse(grepl("R-", V(m)$name), "circle", V(m)$nodeType)

    if (length(interesting_nodes) > 1) {
        V(m)$nodeType <- ifelse(V(m)$name %in% interesting_nodes, "star", V(m)$nodeType)
        #print(V(m)$nodeType)
    }


    m
}


modulePostProcessing <- function(m, reactome, interesting_nodes=NA) {
    reactions <- V(m)[grepl("R-", V(m)$name)]
    vertices <- V(m)[!grepl("R-", V(m)$name)]
    for (i in 1:length(reactions))
        for (j in 1:length(vertices)) {
            r <- reactions[i]$name
            v_name <- vertices[j]$name
            v_label <- vertices[j]$label
            if ((r %in% reactome[symbol==v_name, reaction_id]) && (! r %in% neighbors(m, v_name)$name)) {
                cat("Add edge: ", r, " -- ", v_label, "\n")
                name <- paste(v_label, "is in reactions", r, collapse=" ")
                m <- m + edges(r, v_name, name=name)
            }
        }

    # add shared reactions
    all_edges <- unique(reactome[,.(reaction_id, symbol)])
    v <- V(m)$name[!grepl("R-", V(m)$name)]
    v_pairs <- data.table(t(combn(v, 2)))
    intersect_reactions <- function(V1, V2) intersect(all_edges[symbol==V1, reaction_id],
                                                      all_edges[symbol==V2, reaction_id])
    v_pairs[, common_reactions := mapply(intersect_reactions, V1, V2)]
    additional_reactions <- function(V1, V2, common_reactions) setdiff(common_reactions,
                                                                       intersect(neighbors(m, V1)$name, neighbors(m, V2)$name))
    v_pairs[, add_reactions := mapply(additional_reactions, V1, V2, common_reactions)]

    for (i in 1:nrow(v_pairs)) {
        if (length(unlist(v_pairs[i, add_reactions])) > 0) {
            V1 <- v_pairs[i, V1]
            V2 <- v_pairs[i, V2]
            for (j in 1:length(unlist(v_pairs[i, add_reactions]))) {
                add_r <- unlist(v_pairs[i, add_reactions])[j]
                if (! add_r %in% V(m)$name) {
                    cat("Add vertex: ", add_r, "\n")
                    m <- m + vertex(add_r, label=add_r, comment="Shared reaction from covered Reactome")
                }
                cat("Add edge: ", V1, " -- ", add_r, " -- ", V2, "\n")
                m <- m + edges(c(add_r, V1), name="shared")
                m <- m + edges(c(add_r, V2), name="shared")
            }
        }
    }

    # add reactions url
    V(m)$url <- ifelse(startsWith(V(m)$name, "R-"), sprintf("https://reactome.org/content/detail/%s", V(m)$name), V(m)$url)

    # add labels to edges
    elabel <- unique(reactome[,.(reaction_id, symbol, types_short)])
    elabel <- elabel[, label:=list(typ=paste0(types_short, collapse=";")), by=.(reaction_id, symbol)]
    elabel[, edge:=paste(reaction_id, symbol, sep=" ")]
    elabel[, rev_edge:=paste(symbol, reaction_id, sep=" ")]
    elist <- get.edgelist(m)
    elist <- paste(elist[,1], elist[,2])
    #E(m)$name <- elabel[match(elist, elabel[, edge]), label]
    E(m)$label <- elabel[match(elist, elabel[, edge]), label]
    rev_names <- elabel[match(elist, elabel[, rev_edge]), label] # WTF? Why (gene, reactione) if I add (r, v)
    E(m)$label <- ifelse(!is.na(rev_names), rev_names, E(m)$label)

    # add reaction names
    reaction_names <- data.table(read_csv("inst/reactome/mus_reaction_names.csv"))
    reaction_names[, unique_name:=paste(reaction_name, reaction_id)]
    reaction_names <- reaction_names[match(V(m)$name, reaction_names[, reaction_id]), unique_name]
    V(m)$name <- ifelse(!is.na(reaction_names), reaction_names, V(m)$name)

    # add shape
    is_met <- match(V(m)$name, reactome[, cross_ref])
    V(m)$nodeType <- ifelse(!is.na(is_met), "circle", "square")

    # reactions is circle
    V(m)$nodeType <- ifelse(grepl("R-", V(m)$name), "circle", V(m)$nodeType)

    if (length(interesting_nodes) > 1) {
        V(m)$nodeType <- ifelse(V(m)$name %in% interesting_nodes, "star", V(m)$nodeType)
    }
    m
}



saveModuleToSVG <- function(m, name="M0.vs.M1") {
    dot_file <- paste(name, "dot", sep=".")
    svg_file <- paste(name, "svg", sep=".")
    saveModuleToDot(m, file=dot_file, name=name)
    call <- paste("neato -Tsvg ", dot_file, " > ", svg_file)
    system(call, ignore.stderr = T)
}
