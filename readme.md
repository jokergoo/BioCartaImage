# BioCartaImage


BioCarta is a valuable source of biological pathways which not only provides
well manually curated pathways, but also remarkably intuitive pathway images.
One useful features of pathway analysis which is to highlight genes of
interest on the pathway images has been lost. Since the original source of
BioCarta (biocarte.com) is lost from the internet, we digged out the data from
the internet archive and formatted it into a package.


## Install

```r
devtools::install_github("jokergoo/BioCartaImage")
```

## Usage

```r
library(BioCartaImage)
library(grid)
grid.newpage()
grid.biocarta("h_RELAPathway")
```

![image](https://github.com/jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38/assets/449218/ffc692c2-729f-41cf-a045-bf28168a39c6)


Highlight genes:

```r
grid.newpage()
grid.biocarta("h_RELAPathway", color = c("1387" = "yellow"))
```

![image](https://github.com/jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38/assets/449218/89f3a5d6-3b14-4e6c-a4b4-f7469dbb6528)


Add advanced graphics:

```r
grob = biocartaGrob("h_RELAPathway")

grid.newpage()
grob4 = mark_gene(grob, "1387", function(x, y) {
    pos = pos_by_polygon(x, y)
    pushViewport(viewport(x = pos[1] - 10, y = pos[2], 
        width = unit(4, "cm"), height = unit(4, "cm"), 
        default.units = "native", just = "right"))
    grid.rect(gp = gpar(fill = "red"))
    grid.text("add whatever\nyou want here")
    popViewport()
}, capture = TRUE)
grid.draw(grob4)
```

![image](https://github.com/jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38/assets/449218/cef41abf-ec39-4384-9e51-3e3c809ebac0)



## License

MIT @ Zuguang Gu
