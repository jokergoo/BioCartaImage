#' Pre-computed data objects
#' 
#' @details
#' `BIOCARTA_PATHWAYS`, `PATHWAY2BC`, `PATHWAY2ENTREZ` and `BC2ENTREZ` are collected from
#' web.archive.org (\url{https://web.archive.org/web/20170122225118/https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways}).
#' `PATHWAY2MSIGDB` is collected from MSigDB database (\url{https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=CP:BIOCARTA}).
#' The script for generating these datasets can be found at:
#'
#' ```r
#' system.file("script", "process.R", package = "BioCartaImage")
#' ```
#'
#' @return
#' - `BIOCARTA_PATHWAYS`: A list of pathway objects. The pathway object is explained in [`get_pathway()`].
#' - `PATHWAY2BC`: A two-column data frame of pathway IDs and BC IDs.
#' - `PATHWAY2ENTREZ`: A two-column data frame of pathway IDs and gene Entrez IDs.
#' - `PATHWAY2MSIGDB`: A two-column data frame of pathway IDs and MSigDB IDs.
#' - `BC2ENTREZ`: A two-column data frame of BC IDs and gene EntreZ IDs.
#' 
#' The nodes in the original BioCarta pathways are proteins and some of them do not have one-to-one
#' mapping to genes, such as protein families or complex. Here `BC_ID` is the primary ID of proteins/single nodes
#' in BioCarta Pathways and this package provides mapping to gene EntreZ IDs.
#' 
#' @rdname datasets
"BIOCARTA_PATHWAYS"

#' @rdname datasets
"PATHWAY2BC"

#' @rdname datasets
"PATHWAY2ENTREZ"

#' @rdname datasets
"PATHWAY2MSIGDB"

#' @rdname datasets
"BC2ENTREZ"


.onLoad = function(libname, pkgname) {

    all_vars =  c("BIOCARTA_PATHWAYS", "PATHWAY2BC", "PATHWAY2ENTREZ", "PATHWAY2MSIGDB", "BC2ENTREZ")

    utils::data(list = all_vars, package = pkgname, lib.loc = libname)
    ns = asNamespace(pkgname)

    for(var in all_vars) {
        assign(var, get(var), envir = ns)
        namespaceExport(ns, var)
    }
}
