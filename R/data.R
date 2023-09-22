#' Pre-computated datasets
#' 
#' @details
#' `BIOCARTA_PATHWAYS`: A list of pathways.
#' `PATHWAY2BC`: A two-column data frame.
#' `PATHWAY2ENTREZ`: A two-column data frame.
#' `BC2ENTREZ`: A two-column data frame.
#' 
#' `BCID` is an internal ID type in BioCarta. It is more like an ID for the nodes in the pathways.
#' 
#' @rdname datasets
"BIOCARTA_PATHWAYS"

#' @rdname datasets
"PATHWAY2BC"

#' @rdname datasets
"PATHWAY2ENTREZ"

#' @rdname datasets
"BC2ENTREZ"



.onLoad = function(libname, pkgname) {

    all_vars =  c("BIOCARTA_PATHWAYS", "PATHWAY2BC", "PATHWAY2ENTREZ", "BC2ENTREZ")

    utils::data(list = all_vars, package = pkgname, lib.loc = libname)
    ns = asNamespace(pkgname)

    for(var in all_vars) {
        assign(var, get(var), envir = ns)
        namespaceExport(ns, var)
    }
}
