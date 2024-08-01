#Oncoprint plotting for longitudinal PA
library(ComplexHeatmap)

mat_pa <- read.csv("oncomatrix.csv", row.names = 1)
mat <- as.matrix(mat_pa)

#Colour version
col = c("INDEL" = "indianred", "MUT" = "steelblue3", "FUSION" = "green4", "INV" = "orange2", "SPLICE" = "purple3")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  INDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["INDEL"], col = NA))
  },
  # big red
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  # big green
  FUSION = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["FUSION"], col = NA))
  },
  # big orange
  INV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["INV"], col = NA))
  },
  # big orange
  SPLICE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["SPLICE"], col = NA))
  }
)

column_title = "Grayscale OncoPrint for PA"
heatmap_legend_param = list(title = "Alterations", at = c("INDEL", "FUSION", "MUT", "INV", "SPLICE"), 
                            labels = c("Indel", "Fusion", "SNV", "Inversion", "Splice site"))
names <- colnames(mat)
names <- gsub("_", " ", names)

tiff("PA_Oncoprint_panel.tiff",
     units = "in",
     width = 9,
     height = 5,
     res = 300)
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          pct_side = "right", row_names_side = "left",
          column_order = colnames(mat),
          column_title = NULL, show_column_names = FALSE,
          heatmap_legend_param = heatmap_legend_param,
          row_names_gp = gpar(fontsize = 12),
          bottom_annotation = HeatmapAnnotation(
            text = anno_text(names, rot = 45, offset = unit(1, "npc"), just = "right", gp = gpar(fontsize = 12)),
            annotation_height = max_text_width(colnames(mat))
          )
        )
dev.off()
