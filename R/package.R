
#' The BioCartaImage package
#' 
#' BioCarta is a valuable source of biological pathways which not only provides
#' well manually curated pathways, but also remarkable and intuitive pathway images.
#' One useful features of pathway analysis which is to highlight genes of
#' interest on the pathway images is lost. Since the original source of
#' BioCarta (biocarte.com) is lost from the internet, we digged out the data from
#' the internet archive and formatted it into a package.
#' 
#' The core functionality of this package is to highlight certain genes on the pathway image. 
#' The **BioCartaImage** package wraps the pathway image as well as gene locations 
#' into a graphic object
#' 
#' A simple use is as follows:
#' 
#' ```r
#' library(BioCartaImage)
#' library(grid)
#' grid.newpage()
#' grid.biocarta("h_RELAPathway", color = c("1387" = "yellow"))
#' ```
#' 
#' where `"h_RELAPathway"` is a BioCarta pathway ID, `"1387"` (in the EntreZ ID type) is the gene to be highlighted.
#' [`grid.biocarta()`] is a low-level **grid** graphical function which adds the pathway graphic to a certain
#' position in the plot.
#' 
#' More advanced use is first to create a graphic object (a grob), later to add more complex graphics to it:
#' 
#' ```r
#' grid.newpage()
#' grob = biocartaGrob("h_RELAPathway")
#' grob2 = mark_gene(grob, "1387", function(x, y) {
#'     pos = pos_by_polygon(x, y)
#'     pushViewport(viewport(x = pos[1] - 10, y = pos[2], 
#'         width = unit(4, "cm"), height = unit(4, "cm"), 
#'         default.units = "native", just = "right"))
#'     grid.rect(gp = gpar(fill = "red"))
#'     grid.text("add whatever\nyou want here")
#'     popViewport()
#' }, capture = TRUE)
#' grid.draw(grob2)
#' ```
#' 
#' Here [`biocartaGrob()`] creates a grob for the pathway image and [`mark_gene()`] adds more 
#' graphics which are defined by the self-defined function.
#' 
#' For more details, please go to the vignette of this package.
#' 
#' @name BioCartaImage-package
#' @docType package
NULL
