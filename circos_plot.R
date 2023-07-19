suppressMessages(suppressWarnings(require(ggplot2)))
suppressMessages(suppressWarnings(require(circlize)))
suppressMessages(suppressWarnings(require(grDevices)))
suppressMessages(suppressWarnings(require(ComplexHeatmap)))
suppressMessages(suppressWarnings(require(tidyr)))

scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

netVisual_chord_gene_order_top_label_plot <- function(net = NULL,
                                                 sources.use = NULL, targets.use = NULL,
                                                 lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                                 link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                                 transparency = 0.4, link.border = NA,
                                                 title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                                 thresh = 0.05,
                                                 ...){
  
  df <- net
  df$source <- factor(df$source)
  df$target <- factor(df$target)                                     
  sources.use <- levels(df$source)
  targets.use <- levels(df$target)
  # remove the interactions with zero values
  df <- subset(df, prob > 0)
   #df <- df[order(prob),]
#   df <- df %>%
    # distinct(ligand, receptor, .keep_all = TRUE)
  
#   df <- slice_max(df, n = 30, order_by = prob)
  #df <- slice_max(df, n = 20, order_by = prob)
  
  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }
  
  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  
  cell.order.sources <- levels(df$source)[levels(df$source) %in% sources.use]
  cell.order.targets <- levels(df$target)[levels(df$target) %in% targets.use]
  
  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df <- df[order(df$prob, decreasing = T),]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]
  
  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])
  
  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)
  

  shared <- factor(union(df$source,df$target))
  color.use = scPalette(nlevels(shared))
  names(color.use) <- levels(shared)
  #color.use <- as.character(union(df$source,df$target))
  
  # define edge color
  edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  
  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector
  
  df.plot <- df.ordered.source[,c('ligand', 'receptor', 'prob')]
  
  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  #circos.par(t)
  chordDiagram(df.plot,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight+arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = scale,
               link.target.prop = link.target.prop,
               reduce = reduce,
               ...)
  #circos.info()
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)
  
  
  # https://jokergoo.github.io/circlize_book/book/legends.html
  
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }
  circos.clear()
  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()
  return(gg)
}
