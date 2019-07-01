#setwd("~/bioinf/reactome_network/module_constructor")

makeVertexTable <- function(gene.de.unique, met.de.unique, reactome, eps=0.1, ..., vertex.signal="Sv") {
    v <- gene.de.unique[symbol %in% unique(reactome[,symbol]), .(symbol,score)]

    cat("# data genes in covered reactome", v[,.N], "\n")
    cat("# COMPOUNDs in covered reactome", length(unique(reactome[db_type=="COMPOUND", symbol])), "\n")
    cat("# reactions in covered reactome", length(unique(reactome[, reaction_id])), "\n")

    mets <- unique(reactome[db_type=="COMPOUND", symbol])
    reactions <- unique(reactome[, reaction_id])

    if (!is.null(met.de.unique)) {
        v_met <- met.de.unique[KEGG %in% mets, .(KEGG, score)]
        mets <- setdiff(mets, met.de.unique[KEGG %in% mets, KEGG])
        v <- rbindlist(list(v, v_met))
    }

    v[, signal:=paste("S", rownames(v), sep = "")]

    additional_v <- union(mets, reactions)
    v_eps <- data.table(symbol=additional_v,
                        score=-eps,
                        signal=vertex.signal)
    rbindlist(list(v, v_eps))
}


makeEdgeTable <- function(v.table, reactome, eps=0.1, ..., edge.signal="Se") {
    reactome[, c("signal", "score") := list(edge.signal, -eps)] # TODO by ref :(

    # penalty <- function(types, eps) -eps * length(unlist(strsplit(types, split = ",")))
    # reactome[, score := mapply(penalty, types, eps)]
    # reactome[, signal := paste(edge.signal, rownames(reactome), sep = "")]

    e_table <- reactome[symbol %in% unique(v.table[,symbol]),.(reaction_id, symbol, signal, score)]
    unique(e_table, by=c("reaction_id", "symbol"))
}


getAllMets <- function() {
    metabolites <- keggList("compound")
    metabolites <- data.table(
        metabolite=gsub("cpd:", "", names(metabolites), fixed = T),
        metabolite_name=gsub(";.*$", "", metabolites))
    metabolites[, metabolite_url := sprintf("http://www.kegg.jp/entry/%s", metabolite)]
    setkey(metabolites, metabolite)
    metabolites
}


metabolitesInfo <- function(mets.compound) {
    metabolites <- getAllMets()
    metabolites[metabolite %in% mets.compound,]
}

maskGenes <- function(reactome) { # TODO: grep Ubiquitin in gene names
    # remove Ubiquitin

    #uniprot[, genename := toupper(mapIds(org.Mm.eg.db,
    #                                     mapIds(org.Mm.eg.db, uniprot[,cross_ref], 'ENTREZID', 'UNIPROT'),
    #                                     'GENENAME', 'ENTREZID'))]

    reactome[!grepl("^Ub[a,e,r,c]", symbol), ]
}


maskMets <- function(reactome) {
    mets2mask <- fread(system.file("mets2mask.lst", package="gatom"))$ID
    mets <- getAllMets()
    cation <- mets[grepl("cation", metabolite_name), metabolite]
    strange_guys <- c("C00115", "C00034")
    mets2mask <- Reduce(union, list(mets2mask, cation, strange_guys))
    reactome[!(symbol %in% mets2mask), ]
}


writeGmwcsFiles <- function(gene.de.unique, met.de.unique, reactome, graph.dir="gmwcs", r_eps, e_eps) {
    dir.create(graph.dir, showWarnings = FALSE)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")

    # nodes
    v_table <- makeVertexTable(gene.de.unique, met.de.unique, reactome, eps = r_eps)
    write.table(v_table[,.(symbol, score)], file=nodes.file, sep="\t", row.names=F, quote=F, col.names=F)

    #edges
    e_table <- makeEdgeTable(v_table, reactome, eps = e_eps)
    write.table(e_table[,.(reaction_id, symbol, score)], file=edges.file, sep="\t", row.names=F, quote=F, col.names=F)

    cat("\n# nodes: ", v_table[,.N], "\n")
    cat("# edges: ", e_table[,.N], "\n")

    list(nodes.file=nodes.file, edges.file=edges.file)
}



#' Title
#'
#' @param name
#'
#' @return
#' @import org.Mm.eg.db
#'
#' @examples
load_data <- function(name, replace_with_pval_phosp = F) {
    interesting_nodes <- NA
    if (name == "M0.vs.M1") {
        gene.de.raw <- refine_datatable(read_tsv("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv.gz"))
        gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                       gene.de.raw[,ID],
                                       'SYMBOL', 'REFSEQ')]
        met.de.raw <- refine_datatable(read_tsv("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv.gz"))
    } else if (name == "Th0_Ctrl.vs.Th0_Itac") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/reg_network/gene.de_peg.Th0_Ctrl.vs.Th0_Itac.tsv"))
        gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                       gene.de.raw[, ID],
                                       'SYMBOL', 'ENSEMBL')]
        met.de.raw <- NULL
    } else if (name == "Th17_Ctrl.vs.Th17_Itac") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/reg_network/gene.de_peg.Th17_Ctrl.vs.Th17_Itac.tsv"))
        gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                       gene.de.raw[, ID],
                                       'SYMBOL', 'ENSEMBL')]
        met.de.raw <- NULL
    } else if (name == "WT.vs.NE.ge") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/data/dufour/WT.vs.NE.tsv"))
        met.de.raw <- NULL
    }
    else if (name == "WT.vs.NE.ge_mets") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/data/dufour/WT.vs.NE/WT.vs.NE.tsv"))
        met.de.raw <- refine_datatable(read_tsv("~/bioinf/data/dufour/WT.vs.NE/WT.vs.NE.mets.tsv"))

        p <- data.table(read_tsv("/home/mch/bioinf/data/dufour/phosphoproteome.csv"))
        p <- unique(p[, c("T: Gene names", "N: Student's T-test p-value WT_KIKO")])
        setnames(p, names(p), c("symbol", "pval_phosph"))
        p_unique_pval <- unique(p[order(pval_phosph)], by="symbol")

        if (replace_with_pval_phosp) {
            cat("\nReplace gene.de.raw p_val with min(p_val, p_val_phosph)\n")
            gene.de.raw <- p_unique_pval[gene.de.raw, on="symbol"]
            gene.de.raw[, pval_gene:=pval]
            gene.de.raw[, pval:=pmin(pval, pval_phosph, na.rm = TRUE)]
            print(gene.de.raw)
        }

        print(gene.de.raw)
        interesting_nodes <- p_unique_pval[,symbol]

    }
    else if (name == "WT.vs.KI.ge_mets") {
        gene.de.raw <- read_tsv("~/bioinf/data/dufour/WT.vs.KI/WT.vs.KI.tsv")
        print(gene.de.raw)
        gene.de.raw <- refine_datatable(gene.de.raw)
        met.de.raw <- refine_datatable(read_tsv("~/bioinf/data/dufour/WT.vs.KI/WT.vs.KI.mets.tsv"))

        p <- data.table(read_tsv("/home/mch/bioinf/data/dufour/phosphoproteome.csv"))
        p <- unique(p[, c("T: Gene names", "N: Student's T-test p-value WT_ERBB2KI")])
        setnames(p, names(p), c("symbol", "pval_phosph"))
        p_unique_pval <- unique(p[order(pval_phosph)], by="symbol")

        if (replace_with_pval_phosp) {
            cat("Replace gene.de.raw p_val with min(p_val, p_val_phosph)")
            gene.de.raw <- p_unique_pval[gene.de.raw, on="symbol"]
            gene.de.raw[, pval_gene:=pval]
            gene.de.raw[, pval:=pmin(pval, pval_phosph, na.rm = TRUE)]
            print(gene.de.raw)
        }

        print(gene.de.raw)
        interesting_nodes <- p_unique_pval[,symbol]

    }
    else if (name == "WT.vs.KO.ge_mets") {
        gene.de.raw <- read_tsv("~/bioinf/data/dufour/WT.vs.KO/WT.vs.KO.tsv")
        print(gene.de.raw)
        gene.de.raw <- refine_datatable(gene.de.raw)
        met.de.raw <- refine_datatable(read_tsv("~/bioinf/data/dufour/WT.vs.KO/WT.vs.KO.mets.tsv"))

        p <- data.table(read_tsv("/home/mch/bioinf/data/dufour/phosphoproteome.csv"))
        p <- unique(p[, c("T: Gene names", "N: Student's T-test p-value WT_ERRaKO")])
        setnames(p, names(p), c("symbol", "pval_phosph"))
        p_unique_pval <- unique(p[order(pval_phosph)], by="symbol")

        if (replace_with_pval_phosp) {
            cat("Replace gene.de.raw p_val with min(p_val, p_val_phosph)")
            gene.de.raw <- p_unique_pval[gene.de.raw, on="symbol"]
            gene.de.raw[, pval_gene:=pval]
            gene.de.raw[, pval:=pmin(pval, pval_phosph, na.rm = TRUE)]
            print(gene.de.raw)
        }

        print(gene.de.raw)
        interesting_nodes <- p_unique_pval[,symbol]

    }
    else if (name == "plur.0h.vs.72h.ge") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/data/pluripotency/ge/ge_0h_vs_72h.tsv"))
        gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                                      gene.de.raw[, ID],
                                                      'SYMBOL', 'ENSEMBL')]
        met.de.raw <- NULL
    }
    else if (name == "plur.0h.vs.72h.ge_phosph") {
        gene.de.raw <- refine_datatable(read_tsv("~/bioinf/data/pluripotency/ge/ge_0h_vs_72h.tsv"))
        gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
                                                      gene.de.raw[, ID],
                                                      'SYMBOL', 'ENSEMBL')]
        met.de.raw <- NULL

        p <- data.table(read_tsv("/home/mch/bioinf/data/pluripotency/phosph/phosph_72vs0_limma.tsv"))
        p <- unique(p[, c("ID", "pval")])
        setnames(p, names(p), c("symbol", "pval_phosph"))
        p_unique_pval <- unique(p[order(pval_phosph)], by="symbol")

        cat("\nReplace gene.de.raw p_val with min(p_val, p_val_phosph)\n")
        gene.de.raw <- p_unique_pval[gene.de.raw, on="symbol"]
        gene.de.raw[, pval_gene:=pval]
        gene.de.raw[, pval:=pmin(pval, pval_phosph, na.rm = TRUE)]
        print(gene.de.raw)
        interesting_nodes <- p_unique_pval[,symbol]
    }
    cat("\n# genes in raw data: ", gene.de.raw[,.N], "\n\n")
    list(gene.de.raw=gene.de.raw,
         met.de.raw=met.de.raw,
         interesting_nodes=interesting_nodes)

}


#' process_plur_de_ph
#'
#'
#' @export
process_plur_de_ph <- function() {
    process_dataset <- function(ge_fname, ph_fname, name) {
        cat(ge_fname)
        cat(ph_fname)
        gene.de.raw <- refine_datatable(read_tsv(ge_fname))
        print(head(gene.de.raw))
        setnames(gene.de.raw,
                 names(gene.de.raw),
                 c("symbol", "baseMean", "log2FC", "se", "t", "pval", "pval_adj"))
        # gene.de.raw[, symbol := AnnotationDbi::mapIds(org.Mm.eg.db,
        #                                               gene.de.raw[, ID],
        #                                               'SYMBOL', 'ENSEMBL')]
        met.de.raw <- NULL

        p <- data.table(read_tsv(ph_fname))
        print(head(p))
        setnames(p, names(p),
                 c("symbol", "logFC", "AveExpr", "t", "pval_phosph", "pval_adj", "B"))
        #p <- unique(p[, c("ID", "logFC", "AveExpr", "t", "pval", "pval_adj", "B")])
        #setnames(p, names(p), c("symbol", "pval_phosph"))
        p_unique_pval <- unique(p[order(pval_phosph)], by="symbol")

        cat("\nReplace gene.de.raw p_val with min(p_val, p_val_phosph)\n")
        gene.de.raw <- p_unique_pval[gene.de.raw, on="symbol"]
        gene.de.raw[, pval_gene:=pval]
        gene.de.raw[, pval:=pmin(pval, pval_phosph, na.rm = TRUE)]
        print(gene.de.raw)
        interesting_nodes <- p_unique_pval[,symbol]

        cat("\n# genes in raw data: ", gene.de.raw[,.N], "\n\n")
        data <- list(gene.de.raw=gene.de.raw,
             met.de.raw=met.de.raw,
             interesting_nodes=interesting_nodes)

        name1 <- tail(unlist(strsplit(ge_fname, split = "/")), 1)
        name2 <- tail(unlist(strsplit(ph_fname, split = "/")), 1)
        name <- paste(c(gsub('.{4}$', '', name1), "vs", gsub('.{4}$', '', name2)), collapse = "_")

        run_module_constructor(name=name, data=data)
    }

    ph_res <- list.files(path = "~/bioinf/data/pluripotency/phosph/limma_results", full.names = T)
    ge_res <- list.files(path = "~/bioinf/data/pluripotency/ge/deseq_results", full.names = T)
    modules <- apply(expand.grid(ge_res, ph_res), 1,function(x) process_dataset(x[1], x[2]))
}



#' run_module_constructor
#'
#' @param name
#' @param kpos_gene
#' @param kpos_met
#' @param r
#' @param e
#' @param replace_with_pval_phosp
#' @param data
#'
#' @return
#'
#' @import data.table
#' @import readr
#' @import igraph
#' @import org.Mm.eg.db
#' @import KEGGREST
#' @export

#' @examples
run_module_constructor <- function(name,
                                   kpos_gene = 100,
                                   kpos_met = 100,
                                   r = 0.1,
                                   e = 0.1,
                                   ... ,
                                   replace_with_pval_phosp = F,
                                   data = NA)
{
    if (is.na(data)) {
        data <- load_data(name, replace_with_pval_phosp)
    }
    gene.de.raw <- data$gene.de.raw
    met.de.raw <- data$met.de.raw
    interesting_nodes <- data$interesting_nodes

    print(interesting_nodes)

    gene.de.unique <- prepareData(gene.de.raw, k.positive.scores = kpos_gene)
    cat("# genes after preprocessing: ", gene.de.unique[,.N], "\n\n")

    if (!is.null(met.de.raw)) {
        cat("# mets before preprocessing: ", met.de.raw[,.N], "\n\n")
        met.de.unique <- prepareMet(met.de.raw, kpos_met)
        cat("# mets after preprocessing: ", met.de.unique[,.N], "\n\n")
    } else {
        met.de.unique <- NULL
        cat("Mets data is unavailable \n\n")
    }

    reactome <- data.table(read_csv("inst/reactome/reactome_mus.csv"))
    cat("\n# edges in raw reactome: ", reactome[,.N], "\n")
    cat("# reactions in raw reactome: ", length(unique(reactome[,reaction_id])), "\n\n")
    reactome <- prepareReactome(reactome, gene.de.unique)
    cat("\n# edges in covered reactome: ", reactome[,.N], "\n\n")
    #reactome

    instance <- writeGmwcsFiles(gene.de.unique,
                                met.de.unique,
                                reactome,
                                r_eps = r,
                                e_eps = e)

    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, "/home/mch/bioinf/solver/", sep = ":"))

    cat("\nRunning GMWCS\n\n")
    log.file <- file.path("gmwcs", "log")
    system2("gmwcs",  c("--nodes", instance$nodes.file,
                        "--edges", instance$edges.file,
                        "--threads", 2,
                        "--timelimit", 1200
    ), stdout = log.file, stderr = log.file)


    solution.file <- paste0(instance$edges.file, ".out")
    solution <- data.table(read_tsv(solution.file, col_names = FALSE))
    setnames(solution, names(solution), c("reaction", "symbol", "score"))
    solution <- solution[1:(solution[,.N] - 1),]
    solution <- solution[!grepl("n/a", score),]

    m <- solution2module(solution, gene.de.unique, met.de.unique)
    print(m)
    print(intersect(V(m)$label, interesting_nodes))
    saveModuleToSVG(modulePostProcessing_raw(m, reactome, interesting_nodes),
                    paste0(name, "_kpos_", kpos_gene, "_r_", r, "_e_", e, "_raw"))

    m <- modulePostProcessing(m, reactome, interesting_nodes)
    m_name <- paste0(name, "_kpos_", kpos_gene, "_r_", r, "_e_", e)
    cat("Save module to ", m_name, "\n")
    saveModuleToSVG(m, m_name)
    saveModuleToXgmml(m, paste0(m_name, ".xgmml"))
    m
}



#m_th0 <- run_module_constructor("Th0_Ctrl.vs.Th0_Itac", kpos_gene=100)

#m_th17 <- run_module_constructor("Th17_Ctrl.vs.Th17_Itac")

#m_M <- run_module_constructor("M0.vs.M1")

#WT_vs_NE <- run_module_constructor("WT.vs.NE.ge")

#WT_vs_NE_mets <- run_module_constructor("WT.vs.NE.ge_mets",
#kpos_gene = 200, kpos_met = 200, replace_with_pval_phosp = T)


#m_th0_100 <- run_module_constructor("Th0_Ctrl.vs.Th0_Itac", kpos_gene=100)

#m_th17_100 <- run_module_constructor("Th17_Ctrl.vs.Th17_Itac", kpos_gene=100)

#m_M_100 <- run_module_constructor("M0.vs.M1", kpos_gene = 100, kpos_met = 100)

#WT_vs_NE_100 <- run_module_constructor("WT.vs.NE.ge", kpos_gene = 100)

#WT_vs_NE_mets_100 <- run_module_constructor("WT.vs.NE.ge_mets", kpos_gene = 100)
