
#' All Biocarta pathways
#' 
#' @return A vector of pathway IDs.
#' @export
#' @examples
#' all_pathways()
all_pathways = function() {
	names(BIOCARTA_PATHWAYS)
}

#' Get a single pathway
#' 
#' @param pathway_id A biocarta pathway ID.
#' 
#' @return A `biocarta_pathway` object
#' @export
#' @examples
#' get_pathway("h_RELAPathway")
get_pathway = function(pathway_id) {
	p = BIOCARTA_PATHWAYS[[pathway_id]]
	class(p) = "biocarta_pathway"
	p
}

print.biocarta_pathway = function(x, ...) {
	cat("A BioCarta pathway:\n")
	cat("  ID: ", x$id, "\n", sep = "")
	cat("  Name: ", x$name, "\n", sep = "")
	cat("  ", length(x$genes), " nodes, ", length(unique(BC2ENTREZ$ENTREZ[BC2ENTREZ$BCID %in% x$genes])), " genes", "\n", sep = "")
}

#' Genes in a pathway
#' 
#' @param pathway A pathway ID or a `biocarta_pathway` object.
#' 
#' @return A character vector of Entrez IDs.
#' @export
#' @examples
#' genes_in_pathway("h_RELAPathway")
genes_in_pathway = function(pathway) {
	if(inherits(pathway, "character")) {
		pathway = get_pathway(pathway)
	}
	unique(BC2ENTREZ$ENTREZ[BC2ENTREZ$BCID %in% pathway$genes])
}


entrez_to_BC = function(entrez_id) {
	BC2ENTREZ$BCID[ BC2ENTREZ$ENTREZ %in% entrez_id ]
}

BC_to_entrez = function(bc_id) {
	BC2ENTREZ$ENTREZ[ BC2ENTREZ$BCID %in% bc_id ]
}

.ENV = new.env()
.ENV$image = list()

#' Download the pathway image
#' 
#' @param pathway_id A pathway ID.
#' 
#' @return A `raster` object.
#' @export
#' @importFrom magick image_read
#' @examples
#' img = get_pathway_image("h_RELAPathway")
#' class(img)
get_pathway_image = function(pathway_id) {

	if(inherits(pathway_id, "biocarta_pathway")) {
		pathway_id = pathway_id$id
	}

	if(is.null(.ENV$image[[pathway_id]])) {
		url = paste0("~/project/development/biocarta/image/", pathway_id, ".gif")
		img = image_read(url)
		.ENV$image[[pathway_id]] = as.raster(img)
	}
	.ENV$image[[pathway_id]]
}

#' Draw a BioCarta pathway
#' 
#' @param pathway A pathway ID or a `biocarta_pathway` object.
#' @param color A named vector where names should correspond to Entrez IDs.
#' @param x
#' @param y
#' @param width
#' @param height
#' @param just
#' @param hjust
#' @param vjust
#' @param default.units
#' @param name
#' 
#' @export
#' @rdname pathwayGrob
#' @examples
#' grid.pathway("h_RELAPathway")
grid.pathway = function(pathway, color = NULL, 
	x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
	width = NULL, height = NULL, 
	just = "centre", hjust = NULL, vjust = NULL, 
	default.units = "npc", name = NULL) {

	g = pathwayGrob(pathway = pathway, color = color, x = x, y = y,
		width = width, height = height, just = just, hjust = hjust, vjust = vjust,
		default.units = default.units, name = name)
	grid.draw(g)
}

#' @rdname pathwayGrob
#' @export
#' @import grid
pathwayGrob = function(pathway, color = NULL, 
	x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
	width = NULL, height = NULL, 
	just = "centre", hjust = NULL, vjust = NULL, 
	default.units = "npc", name = NULL) {

	if(inherits(pathway, "character")) {
		pathway = get_pathway(pathway)
	}

	image = get_pathway_image(pathway$id)

	size = dim(image)
	image_height = size[1]
	image_width = size[2]

	shape = pathway$shape
	coords = pathway$coords

	vp = viewport(xscale = c(0, image_width), yscale = c(0, image_height),
		x = x, y = y, default.units = default.units, 
		just = just, hjust = hjust, vjust = vjust, name = name)
	vp$check_size = FALSE
	if(!is.null(width) && is.null(height)) {
		vp$width = width
		vp$height = image_height/image_width * vp$width
	} else if(is.null(width) && !is.null(height)) {
		vp$width = image_width/image_height * height
		vp$height = height
	} else if(!is.null(width) && !is.null(height)) {
		vp$width = width
		vp$height = height
		vp$check_size = TRUE
		vp$original_width = vp$width
		vp$original_height = vp$height
	} else {
		vp$check_size = TRUE
		vp$original_width = vp$width
		vp$original_height = vp$height
	}

	gl = gList(rasterGrob(image))

	n = length(shape)

	gl2 = list()

	genes = pathway$genes
	color2 = NULL
	if(!is.null(color)) {
		for(nm in names(color)) {
			if(nm %in% genes) {
				color2[[nm]] = color[[nm]]
			}
			# if the name is entrez
			bc = BC_to_entrez(nm)
			if(length(bc)) {
				color2[bc] = color[[nm]]
			}
		}
	}
	for(i in seq_len(n)) {
		x = coords[[i]]
		nx = length(x)

		if(shape[i] == "poly") {
			gl2[[i]] = polygonGrob(x[seq_len(nx/2)*2-1], image_height - x[seq_len(nx/2)*2], 
				default.units = "native", gp = gpar(col = color2[[ genes[i] ]]))
		} else if(shape[i] == "rect") {
			gl2[[i]] = rectGrob(x[1], height - x[2], width = x[3] - x[1], height = x[4] - x[2], 
				default.units = "native", just = c("left", "bottom"), gp = gpar(col = color2[[ genes[i] ]]))
		} else if(shape[i] == "circle") {
			gl2[[i]] = circleGrob(x[1], x[2], r = x[3], default.units = "native", gp = gpar(col = color2[[ genes[i] ]]))
		}
	}

	gl = c(gl, gl2)
	class(gl) = "gList"

	gTree(children = gl, vp = vp, cl = "BiocartaPathway")
}


makeContext.BiocartaPathway = function(x) {

	if(x$vp$check_size) {

		width = x$vp$xscale[2]
		height = x$vp$yscale[2]

		vp_w = convertWidth(x$vp$original_width, "in", valueOnly = TRUE)
		vp_h = convertHeight(x$vp$original_height, "in", valueOnly = TRUE)

		if(vp_w/vp_h > width/height) {
			x$vp$width = unit(width/height*vp_h, "in")
			x$vp$height = unit(vp_h, "in")
		} else {
			x$vp$width = unit(vp_w, "in")
			x$vp$height = unit(height/width*vp_w, "in")
		}
	}

	x
}

grobWidth.BiocartaPathway = function(x) {
	if(x$vp$check_size) {

		width = x$vp$xscale[2]
		height = x$vp$yscale[2]

		vp_w = convertWidth(x$vp$original_width, "in", valueOnly = TRUE)
		vp_h = convertHeight(x$vp$original_height, "in", valueOnly = TRUE)

		if(vp_w/vp_h > width/height) {
			x$vp$width = unit(width/height*vp_h, "in")
			x$vp$height = unit(vp_h, "in")
		} else {
			x$vp$width = unit(vp_w, "in")
			x$vp$height = unit(height/width*vp_w, "in")
		}
	}
	x$vp$width
}

grobHeight.BiocartaPathway = function(x) {
	if(x$vp$check_size) {

		width = x$vp$xscale[2]
		height = x$vp$yscale[2]

		vp_w = convertWidth(x$vp$original_width, "in", valueOnly = TRUE)
		vp_h = convertHeight(x$vp$original_height, "in", valueOnly = TRUE)

		if(vp_w/vp_h > width/height) {
			x$vp$width = unit(width/height*vp_h, "in")
			x$vp$height = unit(vp_h, "in")
		} else {
			x$vp$width = unit(vp_w, "in")
			x$vp$height = unit(height/width*vp_w, "in")
		}
	}
	x$vp$height
}
