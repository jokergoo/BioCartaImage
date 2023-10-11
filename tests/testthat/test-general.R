
test_that("test BioCartaImage", {
	expect_error(get_pathway("aaaaa"), "Cannot find pathway")
})

p = get_pathway("h_RELAPathway")
print(p)

get_pathway("h_bArrestinPathway")
get_pathway("bArrestinPathway")
get_pathway("h_barrestinpathway")
get_pathway("barrestinpathway")
get_pathway("BIOCARTA_RELA_PATHWAY")

get_pathway_image("h_bArrestinPathway")
image_dimension("h_bArrestinPathway")

genes = genes_in_pathway("h_bArrestinPathway")
bc_ids = BioCartaImage:::entrez_to_BC(genes)
BioCartaImage:::BC_to_entrez(bc_ids)

grid.newpage()
grid.biocarta("h_RELAPathway")


gb = biocartaGrob("h_RELAPathway", x = 0.2, y = 0.2, just = c("left", "bottom"),
	width = unit(2, "cm"))
grobWidth(gb)
grobHeight(gb)

gb = biocartaGrob("h_RELAPathway", x = 0.2, y = 0.2, just = c("left", "bottom"),
	height = unit(2, "cm"))
grobWidth(gb)
grobHeight(gb)


grid.newpage()
genes = genes_in_pathway("h_RELAPathway")
grid.biocarta("h_RELAPathway", color = c("1387" = "black"))


grid.newpage()
pushViewport(viewport(width = 0.5, height = 0.5))
grid.rect()
grid.biocarta("h_RELAPathway")


grid.newpage()
pushViewport(viewport(width = 0.7, height = 0.1))
grid.biocarta("h_RELAPathway", color = c("1387" = "yellow"))

grid.newpage()
pushViewport(viewport(width = 0.1, height = 0.7))
grid.biocarta("h_RELAPathway", color = c("1387" = "yellow"))

grid.newpage()
grid.biocarta("h_RELAPathway", width = 0.1, height = 0.7, color = c("1387" = "yellow"))
grid.biocarta("h_RELAPathway", height = 0.1, width = 0.7, color = c("1387" = "yellow"))


grid.newpage()
genes = genes_in_pathway("h_RELAPathway")
grob = biocartaGrob("h_RELAPathway")

grob2 = mark_gene(grob, "1387", function(x, y) {
	pos = pos_by_polygon(x, y)
	pointsGrob(pos[1], pos[2], default.units = "native")
})

grob2 = mark_gene(grob, "1387", min_area = 1, function(x, y) {
	pos = pos_by_polygon(x, y)
	pointsGrob(pos[1], pos[2], default.units = "native")
})

grob3 = mark_gene(grob, "1387", function(x, y) {
	pos = pos_by_polygon(x, y)
	grid.points(pos[1], pos[2], default.units = "native")
}, capture = TRUE)

grob4 = mark_gene(grob, "1387", function(x, y) {
	pos = pos_by_polygon(x, y)
	pushViewport(viewport(x = pos[1], y = pos[2], width = unit(4, "cm"), height = unit(4, "cm"), 
		default.units = "native", just = "right"))
	grid.rect()
	popViewport()
}, capture = TRUE)
