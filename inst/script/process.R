
setwd("~/project/development/BioCartaImage")

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

	html = read_html(link)

	image_link = html %>% html_nodes(xpath = "//img[@name]") %>% html_attr("src")

	area = html %>% html_nodes(xpath = "//img[@name]/following-sibling::map/area")
	shape = area %>% html_attr("shape")
	coords = area %>% html_attr("coords")
	bc_link = area %>% html_attr("href")

	list(
		image_link = image_link,
		shape = shape,
		coords = coords,
		bc_link = bc_link
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
	l = grepl("/Genes/", x$bc_link)
	x$shape = x$shape[l]
	x$coords = x$coords[l]
	x$coords = lapply(strsplit(x$coords, ","), as.numeric)
	x$bc_link = x$bc_link[l]
	x$bc_link = gsub(" ", "%20", x$bc_link)
	x$bc = gsub("^.*BCID=(.*)$", "\\1", x$bc_link)
	x
})

pathway_list = lapply(pl, function(x) {
	x2 = x[c("id", "name", "bc", "shape", "coords")]
	x2$image_file = basename(x$image_link)
	class(x2) = "biocarta_pathway"
	x2
})
names(pathway_list) = sapply(pathway_list, function(x) x$id)



BIOCARTA_PATHWAYS = pathway_list

html = read_html("https://data.broadinstitute.org/gsea-msigdb/msigdb/biocarta/human/")
msigdb_biocarta_list = html %>% html_element("table") %>% html_table()
msigdb_biocarta_list = msigdb_biocarta_list$Name
msigdb_biocarta_list = msigdb_biocarta_list[grep("gif$", msigdb_biocarta_list)]

for(nm in names(pathway_list)) {
	i = which(tolower(msigdb_biocarta_list) == tolower(pathway_list[[nm]]$image_file))
	if(length(i)) {
		pathway_list[[nm]]$msigdb_image_file = msigdb_biocarta_list[i]
	} else {
		pathway_list[[nm]]$msigdb_image_file = ""
	}
}
save(BIOCARTA_PATHWAYS, file = "data/BIOCARTA_PATHWAYS.RData", compress = "xz")

# download all gif images
for(i in seq_along(pl)) {
	download.file(pl[[i]]$image_link, destfile = basename(pl[[i]]$image_link))
}

all_bc_links = unique(unlist(lapply(pl, function(x) x$bc_link)))
all_bc_links = all_bc_links[grepl("/Genes/", all_bc_links)]
all_bc_links = gsub(" ", "%20", all_bc_links)

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
for(i in i:length(all_bc_links)) {
	gl[[i]] = read_gene(all_bc_links[i], prefix = i)
}

# validate failed links
for(i in which(sapply(gl, function(x) identical(x, "")))) {
	gl[[i]] = read_gene(all_bc_links[i], prefix = i)
}

gl = lapply(gl, function(x) unique(x[x != ""]))
gl = tapply(seq_along(gl), gsub("^.*BCID=(.*)$", "\\1", all_bc_links), function(ind) {
	as.character(unique(unlist(gl[ind])))
})

BC2ENTREZ = data.frame(BCID = rep(names(gl), times = sapply(gl, length)),
	                   ENTREZ = unlist(gl))

lt_foo = lapply(pathway_list, function(x) unique(x$bc))
PATHWAY2BC = data.frame(PATHWAY = rep(sapply(pathway_list, function(x) x$id), times = sapply(lt_foo, length)),
	                    BCID = unlist(lt_foo))

lt_foo2 = lapply(pathway_list, function(x) unique(unlist(gl[x$bc])))
PATHWAY2ENTREZ = data.frame(PATHWAY = rep(sapply(pathway_list, function(x) x$id), times = sapply(lt_foo2, length)),
	                    ENTREZ = unlist(lt_foo2))



## get mapping betweeo MsigDb IDs and biocarta pathway IDs
html = read_html("https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=CP:BIOCARTA")
msigdb_all_pathways = html %>% html_elements("#geneSetTable td a") %>% html_attr("href")
msigdb2biocarta = NULL
for(i in seq_along(msigdb_all_pathways)) {
	cat(msigdb_all_pathways[[i]], "\n")
	html = read_html(paste0("https://www.gsea-msigdb.org/gsea/", msigdb_all_pathways[i]))
	tb = html %>% html_element(".lists4") %>% html_table()
	tb[tb[, 1] == "External links", 2]
	

	msigdb_id = gsub("\\.html$", "", basename(msigdb_all_pathways[i]))
	biocarta_id = gsub("\\.gif$", "", basename(tb[[2]][tb[[1]] == "External links"]))

	msigdb2biocarta[[biocarta_id]] = msigdb_id
}
PATHWAY2MSIGDB = data.frame(MSIGDB = unname(unlist(msigdb2biocarta)), BIOCARTA = names(msigdb2biocarta))

save(BC2ENTREZ, file = "BC2ENTREZ.RData", compress = "xz")
save(PATHWAY2BC, file = "PATHWAY2BC.RData", compress = "xz")
save(PATHWAY2ENTREZ, file = "PATHWAY2ENTREZ.RData", compress = "xz")
save(PATHWAY2MSIGDB, file = "PATHWAY2MSIGDB.RData", compress = "xz")
