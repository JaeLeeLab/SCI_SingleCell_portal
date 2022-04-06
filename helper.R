# Helper functions (e.g. plotting) for the website

load(file = './data/appdata.RData')

f <- h5file(filename = './data/sci.h5',
            mode = 'r')
log_x_sci <- vector(mode = 'list', length = 4000)
log_x_sci_timestamp <- vector(mode = 'list', length = 4000)
data_index <- 1


# Helper functions --------------------------------------------------------

load_data <- function(
  feature
) {
  if (!feature %in% names(log_x_sci)) {
    if (!all(sapply(log_x_sci, is.null))) {
      data_index <<- which.min(log_x_sci_timestamp)
    }
    g.i <- which(f[[names(f)]][['barcodes']][1:nrow(vars_sci)] == feature)
    i <- c(f[[names(f)]][['indptr']][g.i:(g.i + 1)]) + 1
    nonzero_vals <- f[[names(f)]][['data']][i[1]:(i[2]-1)]
    c <- f[[names(f)]][['indices']][i[1]:(i[2]-1)] + 1
    tmp <- rep.int(0, times = ncells)
    tmp[c] <- nonzero_vals
    tmp <- log1p(tmp/library_size * 10000)
    log_x_sci[[data_index]] <<- tmp
    names(log_x_sci)[data_index] <<- feature
    log_x_sci_timestamp[[data_index]] <<- Sys.time()
    names(log_x_sci_timestamp[[data_index]]) <<- feature
    data_index <<- which(sapply(log_x_sci, is.null))[1]
  }
}

get_feature_colorscale <- function(x) {
  a <- x/max(x)
  cols <- rgb(colorRamp(colors = c('grey90','skyblue','royalblue','darkorchid','violetred','red'), bias = 0.9)(a), maxColorValue = 256)
  return(cols)
}

# shuffle <- function(x) {
#   return(x[sample(1:length(x), length(x))])
# }
# 
# partial_order <- function(x, how.much = 1, decreasing = FALSE) {
#   degree <- ceiling(length(x)/how.much)
#   x <- split(order(x, decreasing = decreasing), rep(1:degree, each = how.much)[1:length(x)])
#   x <- unlist(lapply(x, shuffle), use.names = FALSE, recursive = FALSE)
#   return(x)
# }
# 
# partial_order_positive <- function(x, how.much, decreasing = FALSE) {
#   pos <- x[x > 0]
#   not <- which(x <= 0)
#   if (decreasing) {
#     return(c(not, partial_order(pos, how.much = how.much, decreasing = decreasing)))
#   } else {
#     return(c(partial_order(pos, how.much = how.much, decreasing = decreasing), not))
#   }
# }

draw_dimplot <- function(
  dataset_value, 
  groupby,
  ptsize,
  xranges,
  yranges,
  draw_labels
) {
  plot_these <- c(paste0(dataset_value, c('_UMAP_1', '_UMAP_2')), groupby)
  umap_df <- obs_sci[plot_these]
  colnames(umap_df) <- c('dim1', 'dim2', 'groupby')
  umap_df <- umap_df[!is.na(umap_df$dim1),]
  if (is.factor(umap_df$groupby)) {
    umap_df$groupby <- droplevels(umap_df$groupby)
  }
  umap_df <- umap_df[sample(1:nrow(umap_df), nrow(umap_df)),]
  label_df <- label_coords[plot_these]
  colnames(label_df) <- c('dim1', 'dim2', 'label')
  label_df <- label_df[!is.na(label_df$dim1) & !is.na(label_df$label),]
  label_df$label <- factor(label_df$label, levels = levels(umap_df$groupby))
  n <- nrow(label_df)
  label_cex <- (n^2 + n + 1) / (n^2 - n + 1)
  
  t1 <- Sys.time()
  ggplot(data = umap_df) + 
    geom_point(mapping = aes(x = dim1, y = dim2, color = groupby),
               size = 0.5)
  t1 <- Sys.time() - t1
  t2 <- Sys.time()
  plot(x = umap_df$dim1, y = umap_df$dim2, 
       col = rainbow(length(unique(umap_df$groupby)), v = 0.9)[umap_df$groupby],
       pch = 16, cex = ptsize, xlab = 'UMAP_1', ylab = 'UMAP_2',
       main = paste('Grouped by:', groupby),
       xaxt = 'n', yaxt = 'n', cex.main = 1.5,
       xlim = xranges, ylim = yranges
  )
  t2 <- Sys.time() - t2  
  
  {
    par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2,0.2,2,0.2), xpd = FALSE)
    plot(x = umap_df$dim1, y = umap_df$dim2, 
         col = rainbow(length(unique(umap_df$groupby)), v = 0.9)[umap_df$groupby],
         pch = 16, cex = ptsize, xlab = 'UMAP_1', ylab = 'UMAP_2',
         main = paste('Grouped by:', groupby),
         xaxt = 'n', yaxt = 'n', cex.main = 1.5,
         xlim = xranges, ylim = yranges
    )
    if (draw_labels) {
      text(x = label_df$dim1,
           y = label_df$dim2,
           labels = label_df$label,
           cex = label_cex)
    }
  }
}

draw_featureplot <- function(
  dataset_value, 
  feature,
  ptsize,
  xranges,
  yranges
) {
  umap_df <- obs_sci[paste0(dataset_value, c('_UMAP_1', '_UMAP_2'))]
  load_data(feature = feature)
  umap_df$feature <- log_x_sci[[feature]]
  colnames(umap_df) <- c('dim1','dim2','feature')
  umap_df <- umap_df[order(umap_df$feature),]
  par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2,0.2,2,0.2), xpd = FALSE)
  if (!sum(umap_df$feature) > 0) {
    plot(x = umap_df$dim1, y = umap_df$dim2,
         col = 'grey90',
         pch = 16, cex = ptsize, xlab = 'UMAP_1', ylab = 'UMAP_2',
         xaxt = 'n', yaxt = 'n', cex.main = 1.5,
         main = paste('Gene:', feature),
         xlim = xranges, ylim = yranges)
  } else {
    plot(x = umap_df$dim1, y = umap_df$dim2,
         col = get_feature_colorscale(umap_df$feature),
         pch = 16, cex = ptsize, xlab = 'UMAP_1', ylab = 'UMAP_2',
         xaxt = 'n', yaxt = 'n', cex.main = 1.5,
         main = paste('Gene:', feature),
         xlim = xranges, ylim = yranges)
  }
}

splitsubplot <- function(
  umap_df, 
  time, 
  ptsize,
  xranges,
  yranges
) {
  tmp_df <- umap_df[umap_df$time == time,]
  plot(x = tmp_df$dim1, y = tmp_df$dim2,
       col = tmp_df$col, pch = 16, cex = ptsize*1.5, 
       xlab = 'UMAP_1', ylab = 'UMAP_2',
       xaxt = 'n', yaxt = 'n', main = paste('Time:', time), cex.main = 2,
       xlim = xranges, ylim = yranges)
}

draw_splitfeatureplot <- function(
  dataset_value, 
  feature,
  ptsize,
  xranges,
  yranges
) {
  umap_df <- obs_sci[c(paste0(dataset_value, c('_UMAP_1', '_UMAP_2')), 'InjuryTimePoint')]
  load_data(feature = feature)
  umap_df$feature <- log_x_sci[[feature]]
  colnames(umap_df) <- c('dim1','dim2','time', 'feature')
  if (!sum(umap_df$feature) > 0) {
    umap_df$col <- 'grey90'
  } else {
    umap_df$col <- get_feature_colorscale(umap_df$feature)
  }
  umap_df <- umap_df[order(umap_df$feature),]
  par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2,0.2,2,0.2), xpd = FALSE, mfrow = c(1,4))
  {
    splitsubplot(umap_df = umap_df, time = levels(umap_df$time)[1], ptsize,
                 xranges = xranges, yranges = yranges)
    splitsubplot(umap_df = umap_df, time = levels(umap_df$time)[2], ptsize,
                 xranges = xranges, yranges = yranges)
    splitsubplot(umap_df = umap_df, time = levels(umap_df$time)[3], ptsize,
                 xranges = xranges, yranges = yranges)
    splitsubplot(umap_df = umap_df, time = levels(umap_df$time)[4], ptsize,
                 xranges = xranges, yranges = yranges)
  }
}

draw_dimplotlegend <- function(
  dataset_value, 
  groupby,
  ptsize,
  labelsize
) {
  plot_these <- c(paste0(dataset_value, c('_UMAP_1', '_UMAP_2')), groupby)
  label_df <- label_coords[plot_these]
  colnames(label_df) <- c('dim1','dim2','label')
  label_df$label <- factor(label_df$label, levels = levels(obs_sci[[groupby]]))
  label_df <- label_df[!is.na(label_df$dim1) & !is.na(label_df$label),]
  label_df$label <- droplevels(label_df$label)
  legend_labels <- levels(label_df$label)
  n <- length(legend_labels)
  label_y_spacing <- (n^2 + n + 5) / (n^2 - n + 5) - 0.060337
  par(mai = c(0.1, 0.1, 1, 0.1), mar = c(0.2, 0.2, 2, 0.2), xpd = FALSE)
  {
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend(x = 'topleft',
           y.intersp = 1.05,
           legend = legend_labels,
           col = rainbow(length(legend_labels), v = 0.9)[1:length(legend_labels)],
           pch = 16, xpd = TRUE, ncol = 1, plot = TRUE,
           pt.cex = ptsize*0.5, 
           # cex = labelsize, 
           bty = 'n')
  }
}

draw_featuredotplot <- function(
  dataset_value,
  feature,
  groupby
) {
  plot_these <- c(paste0(dataset_value, '_UMAP_1'), groupby, 'InjuryTimePoint')
  data_df <- obs_sci[plot_these]
  load_data(feature = feature)
  data_df$feature <- log_x_sci[[feature]]
  if (!sum(data_df$feature) > 0) {
    text(x = 0.5, y = 0.5, 'Gene not detected in any cells.',
         cex = 1.5, col = 'black')
    return()
  }
  colnames(data_df) <- c('dataset', 'groupby', 'time', 'feature')
  data_df <- data_df[!is.na(data_df$groupby) & 
                       !is.na(data_df$feature) &
                       !is.na(data_df$dataset),]
  pct <- sapply(
    X = split(
      x = data_df$feature, 
      f = paste(data_df$groupby, data_df$time, sep = '.')
    ),
    FUN = function(x) sum(x > 0)/length(x) * 100
  )
  avg <- sapply(
    X = split(
      x = data_df$feature, 
      f = paste(data_df$groupby, data_df$time, sep = '.')
    ),
    FUN = mean
  )
  tmp_group <- sapply(strsplit(names(pct), '\\.'), `[`, 1)
  tmp_group <- factor(
    x = tmp_group, 
    levels = rev(levels(droplevels(obs_sci[[groupby]])))
  )
  tmp_time <- sapply(strsplit(names(pct), '\\.'), `[`, 2)
  tmp_ncells <- as.numeric(table(paste(data_df$groupby, data_df$time, sep = '.')))
  feature_table <- data.frame(
    'ncells' = tmp_ncells,
    'pct' = pct,
    'avg' = avg,
    'groupby' = tmp_group,
    'time' = tmp_time
  )
  feature_table$time <- factor(
    x = feature_table$time,
    levels = levels(obs_sci$InjuryTimePoint)
  )
  feature_limit <- ceiling(max(avg) * 10) / 10
  
  scale_text <- (length(unique(feature_table$groupby)) + 5)/length(unique(feature_table$groupby))
  ggplot(data = feature_table, 
         mapping = aes(x = time, y = groupby)) +
    geom_point(mapping = aes(size = pct, fill = avg), pch = 21)+
    scale_fill_viridis_c(option = 'A',
                         breaks = c(0, feature_limit),
                         limits = c(0, feature_limit),
                         labels = c(0, feature_limit)) +
    scale_size(range = c(0, 8), limits = c(0, 100)) +
    theme_bw() +
    xlab(label = 'Time after injury') +
    theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                      size = 16),
          legend.title.align = 1,
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.spacing.y = unit(1, "cm")) +
    guides(fill = guide_colorbar(title = 'log-normalized\ncounts',
                                 title.position = 'left',
                                 frame.colour = 'black',
                                 title.hjust = 0.5,
                                 title.vjust = 0.5,
                                 ticks = FALSE),
           size = guide_legend(title = 'Percent detected',
                               title.position = 'left',
                               override.aes = list(color = 'black',
                                                   fill = 'black')))
}