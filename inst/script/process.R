


library(rvest)

# biocarta.com has a copy on cgap.nci.nih.gov, which is archived on archive.org

html = read_html("https://web.archive.org/web/20170122225118/https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways")
al = html %>% html_elements("#content table li a.genesrch")

links = al %>% html_attr("href")
names = al %>% html_text()
l = names != ""
links = links[l]
names = names[l]

read_pathway = function(link) {
	cat(link, "\n")
	Sys.sleep(1)

	html = read_html(link)

	image_link = html %>% html_nodes(xpath = "//img[@name]") %>% html_attr("src")

	area = html %>% html_nodes(xpath = "//img[@name]/following-sibling::map/area")
	shape = area %>% html_attr("shape")
	coords = area %>% html_attr("coords")
	gene_link = area %>% html_attr("href")

	list(
		image_link = image_link,
		shape = shape,
		coords = coords,
		gene_link = gene_link
	)
}

# a list of pathways, images and links to genes
pl = list()
for(i in 1:length(links)) {
	pl[[i]] = read_pathway(links[[i]])
	pl[[i]]$link = links[i]
	pl[[i]]$name = names[i]
	pl[[i]]$id = basename(links[i])
}

# BC id, or node id
pl = lapply(pl, function(x) {
	l = grepl("/Genes/", x$gene_link)
	x$shape = x$shape[l]
	x$coords = x$coords[l]
	x$coords = lapply(strsplit(x$coords, ","), as.numeric)
	x$gene_link = x$gene_link[l]
	x$gene_link = gsub(" ", "%20", x$gene_link)
	x$genes = gsub("^.*BCID=(.*)$", "\\1", x$gene_link)
	x
})

pathway_list = lapply(pl, function(x) {
	x[c("id", "name", "genes", "shape", "coords")]
})
names(pathway_list) = sapply(pathway_list, function(x) x$id)

BIOCARTA_PATHWAYS = pathway_list
save(BIOCARTA_PATHWAYS, file = "BIOCARTA_PATHWAYSRData", compress = "xz")

# download all gif images
for(i in seq_along(pl)) {
	download.file(pl[[i]]$image_link, destfile = basename(pl[[i]]$image_link))
}

all_gene_links = unique(unlist(lapply(pl, function(x) x$gene_link)))
all_gene_links = all_gene_links[grepl("/Genes/", all_gene_links)]
all_gene_links = gsub(" ", "%20", all_gene_links)

read_gene = function(link, prefix = "") {
	cat(prefix, link, "\n")
	oe = try(html <- read_html(link))
	if(inherits(oe, "try-error")) {
		if(grepl("LLNO=", link)) {
			return(gsub("^.*LLNO=(\\d+).*$", "\\1", link))
		}
	}

	type = html %>% html_elements("#content h3") %>% html_text() 
	if(length(type) == 0) {
		if(grepl("LLNO=", link)) {
			return(gsub("^.*LLNO=(\\d+).*$", "\\1", link))
		} else {
			return("")
		}
	}
	if(type == "Gene List") {
		gl = html %>% html_nodes(xpath = "//div[@id='content']/table[last()]/tr/td[last()]/a")
		if(length(gl) == 0) {
			tb = html %>% html_element("#content table") %>% html_table()
			l = tb[, 1] == "GeneFinder Results For:"
			if(any(l)) {
				return(unique(strsplit(gsub("^.*;\\s*", "", tb[l, 2]), ",")[[1]]))
			} else {
				return("")
			}
		}
		g = gl %>% html_attr("href")
		g = unique(g)
		g = paste0(dirname(link), "/", g)
		return(unname(unlist(lapply(g, read_gene, prefix = "  "))))
	}

	tb = html %>% html_element("#content table") %>% html_table()
	if(length(tb) == 0) {
		return("")
	}
	tb[, 2][ tb[,1] == "Entrez Gene ID:" ]
}

## genes: biocarta ID to entrez ID
# gl = list()
for(i in i:length(all_gene_links)) {
	gl[[i]] = read_gene(all_gene_links[i], prefix = i)
}

# validate failed links
for(i in which(sapply(gl, function(x) identical(x, "")))) {
	gl[[i]] = read_gene(all_gene_links[i], prefix = i)
}

gl = lapply(gl, function(x) unique(x[x != ""]))
gl = tapply(seq_along(gl), gsub("^.*BCID=(.*)$", "\\1", all_gene_links), function(ind) {
	as.character(unique(unlist(gl[ind])))
})

BC2ENTREZ = data.frame(BCID = rep(names(gl), times = sapply(gl, length)),
	                   ENTREZ = unlist(gl))

lt_foo = lapply(pathway_list, function(x) unique(x$genes))
PATHWAY2BC = data.frame(PATHWAY = rep(sapply(pathway_list, function(x) x$id), times = sapply(lt_foo, length)),
	                    BCID = unlist(lt_foo))

lt_foo2 = lapply(pathway_list, function(x) unique(unlist(gl[x$genes])))
PATHWAY2ENTREZ = data.frame(PATHWAY = rep(sapply(pathway_list, function(x) x$id), times = sapply(lt_foo2, length)),
	                    ENTREZ = unlist(lt_foo2))

save(BC2ENTREZ, file = "BC2ENTREZ.RData", compress = "xz")
save(PATHWAY2BC, file = "PATHWAY2BC.RData", compress = "xz")
save(PATHWAY2ENTREZ, file = "PATHWAY2ENTREZ.RData", compress = "xz")
