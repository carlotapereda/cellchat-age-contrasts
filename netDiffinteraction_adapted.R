netVisual_diffInteraction_adapted <- function(object,
                                       comparison = c(1, 2),
                                       measure = c("count", "weight", "count.merged", "weight.merged"),
                                       color.use = NULL,
                                       color.edge = c("#b2182b", "#2166ac"),
                                       title.name = NULL,
                                       sources.use = NULL,
                                       targets.use = NULL,
                                       remove.isolate = FALSE,
                                       top = 1,
                                       weight.scale = FALSE,
                                       vertex.weight = 20,
                                       vertex.weight.max = NULL,
                                       vertex.size.max = 15,
                                       vertex.label.cex = 1,
                                       vertex.label.color = "black",
                                       edge.weight.max = NULL,
                                       edge.width.max = 8,
                                       alpha.edge = 0.6,
                                       label.edge = FALSE,
                                       edge.label.color = "black",
                                       edge.label.cex = 0.8,
                                       edge.curved = 0.2,
                                       shape = "circle",
                                       layout = in_circle(),
                                       margin = 0.2,
                                       arrow.width = 1,
                                       arrow.size = 0.2)


{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, 
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  #igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$color <- "grey"
  igraph::V(g)$frame.color <- "grey"
  #igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  
  # Abbreviation map
  CELLTYPE_MAP <- c(
    "Astrocytes"                = "Ast",
    "Excitatory.Neurons"        = "Ex.Neu",
    "Inhibitory.Neurons"        = "In.Neu",
    "Microglia"                 = "Mic",
    "Oligodendrocyte.Precursor" = "OPCs",
    "Oligodendrocytes"          = "Oli"
  )
  
  # Old:
  # vertex_labels <- c("Astrocytes", "Excitatory.Neurons", "Inhibitory.Neurons", "Microglia", "Oligodendrocyte.Precursor", "Oligodendrocytes")
  # V(g)$label <- vertex_labels
  
  # New:
  # Match vertex names to CELLTYPE_MAP (if present, otherwise keep original)
  vertex_labels <- ifelse(V(g)$name %in% names(CELLTYPE_MAP),
                          CELLTYPE_MAP[V(g)$name],
                          V(g)$name)
  
  V(g)$label <- vertex_labels
  
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}
