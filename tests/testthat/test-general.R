
test_that("test BioCartaImage", {
	expect_error(get_pathway("aaaaa"), "Cannot find pathway")
})


get_pathway("h_bArrestinPathway")
get_pathway("bArrestinPathway")
get_pathway("h_barrestinpathway")
get_pathway("barrestinpathway")

get_pathway_image("h_bArrestinPathway")

grid.newpage()
grid.biocarta("h_RELAPathway")


grid.newpage()
grid.biocarta("h_RELAPathway", x = 0.2, y = 0.2, just = c("left", "bottom"),
	width = unit(6, "cm"))


grid.newpage()
genes = genes_in_pathway("h_RELAPathway")
grid.biocarta("h_RELAPathway", color = c("1387" = "black"))


grid.newpage()
pushViewport(viewport(width = 0.5, height = 0.5))
grid.rect()
grid.biocarta("h_RELAPathway")




grid.newpage()
genes = genes_in_pathway("h_RELAPathway")
grob = biocartaGrob("h_RELAPathway")

grob2 = mark_gene(grob, "1387", function(x, y) {
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
