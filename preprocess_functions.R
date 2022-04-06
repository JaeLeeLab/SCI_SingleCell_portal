

# setup packages ----------------------------------------------------------

setup_packages <- function(
  required.packages = c('Seurat','dplyr','HDF5Array', 'Matrix', 'plyr', 'tidyr'),
  additional.packages = NULL
) {
  all.packages <- c(required.packages, additional.packages)
  load.successful <- sapply(
    X = all.packages, 
    FUN = require, 
    character.only = TRUE
  )
  if (any(!load.successful)) {
    needs.install <- all.packages[!load.successful]
    error.message <- paste('the following packages need to be installed:', needs.install, sep = ' ')
    stop(error.message)
  }
}


# util function ------------------------------------------------------------

factors_to_characters <- function(df) {
  rnames <- rownames(df)
  converted <- data.frame(lapply(
    X = df, 
    FUN = function(x) if(is.factor(x)) as.character(x) else x
  ))
  rownames(converted) <- rnames
  return(converted)
}

# extract cell metadata ----------------------------------------------------

extract_cell_metadata <- function(
  seurat.object,
  reduction = 'umap'
) {
  library_size <- seurat.object@meta.data$nCount_RNA
  ncells <- ncol(seurat.object)
  coldata <- factors_to_characters(seurat.object@meta.data)
  
  dimreduc <- slot(seurat.object[[reduction]], 'cell.embeddings')
  colnames(dimreduc) <- paste(reduction, c('1', '2'), sep = '_')
  return(list('library_size' = library_size,
              'ncells' = ncells, 
              'coldata' = coldata, 
              'dimreduc' = dimreduc))
}


# extract feature metadata ------------------------------------------------

extract_feature_metadata <- function(
  seurat.object,
  default.assay
) {
  rowdata <- slot(seurat.object[[default.assay]], 'meta.features')
  # potentially include stuff for variable features
  return(list('rowdata' = rowdata))
}


# Extract seurat data -----------------------------------------------------

# Extract counts matrix, cell (column) meta-data, and feature (row) meta-data from a given seurat object.
extract_seurat_data <- function(
  seurat.object.path,
  include.counts = TRUE,
  default.assay = 'RNA',
  reduction = 'umap'
) {
  message('Processing: ', seurat.object.path)
  message('reading in...')
  seurat.object <- readRDS(file = seurat.object.path)
  if (class(seurat.object) != 'Seurat') {
    stop('seurat.object is not a Seurat Object')
  }
  # Discard everything but metadata and default.assay count matrices
  DefaultAssay(seurat.object) <- default.assay
  seurat.object <- DietSeurat(
    object = seurat.object,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL, # all features
    assays = default.assay,
    dimreducs = reduction,
    graphs = NULL
  )
  message('extracting metadata and counts...')
  seurat.object.cell.metadata <- extract_cell_metadata(
    seurat.object = seurat.object,
    reduction = reduction
  )
  seurat.object.feature.metadata <- extract_feature_metadata(
    seurat.object = seurat.object,
    default.assay = default.assay
  )
  outs <- list('cell.metadata' = seurat.object.cell.metadata,
               'feature.metadata' = seurat.object.feature.metadata)
  if (include.counts) {
    outs[['counts']] <- slot(seurat.object[[default.assay]], 'counts')
  }
  return(outs)
}

# aggregate cell coldata --------------------------------------------------

aggregate_metadata <- function(
  full.data,
  subset.data
) {
  subset.coldata <- lapply(
    X = subset.data,
    FUN = function(x) {
      tmp <- x[['cell.metadata']][['coldata']]
      tmp$rn <- rownames(tmp)
      return(tmp)
    }
  )
  combined.coldata <- c(
    list(full.data[['cell.metadata']][['coldata']]),
    subset.coldata
  )
  shared.cols <- Reduce(
    f = intersect,
    x = lapply(X = combined.coldata, FUN = colnames)
  )
  subset.coldata <- Reduce(
    f = function(x, y) dplyr::full_join(x, y),
    x = subset.coldata
  )
  full.data[['cell.metadata']][['coldata']]$rn <- rownames(
    full.data[['cell.metadata']][['coldata']]
  )
  combined.coldata <- dplyr::full_join(
    x = full.data[['cell.metadata']][['coldata']],
    y = subset.coldata,
    by = c('rn')
  )
  combined.coldata <- combined.coldata[!grepl(
    pattern = '\\.y',
    x = colnames(combined.coldata)
  )]
  colnames(combined.coldata) <- gsub('\\.x', '', colnames(combined.coldata))
  # rownames(combined.coldata) <- combined.coldata$rn
  # combined.coldata$rn <- NULL
  return(combined.coldata)
}


# aggregate dimreduction coordinates --------------------------------------

aggregate_dimreduc <- function(
  full.data,
  subset.data,
  full.data.prefix,
  subset.data.prefixes
) {
  combined.dimreduc <- c(
    list('full' = as.data.frame(full.data[['cell.metadata']][['dimreduc']])),
    lapply(
      X = subset.data,
      FUN = function(x) as.data.frame(x[['cell.metadata']][['dimreduc']])
    )
  )
  combined.prefixes <- c(full.data.prefix, subset.data.prefixes)
  if (length(combined.prefixes) != length(combined.dimreduc)) {
    stop('length of prefixes does not match length of subset.data')
  }
  for (i in 1:length(combined.dimreduc)) {
    colnames(combined.dimreduc[[i]]) <- paste(
      combined.prefixes[i],
      colnames(combined.dimreduc[[i]]),
      sep = '_'
    )
    combined.dimreduc[[i]]$rn <- rownames(combined.dimreduc[[i]])
  }
  combined.dimreduc <- Reduce(
    f = function(x, y) dplyr::full_join(x, y, by = 'rn'),
    x = combined.dimreduc
  )
  # rownames(combined.dimreduc) <- combined.dimreduc$rn
  # combined.dimreduc$rn <- NULL
  return(combined.dimreduc)
}

preprocess_seurat_data <- function(
  full.data.path,
  full.data.prefix = 'NewStudy',
  subset.data.paths = NULL,
  subset.data.prefixes = NULL,
  default.assay = 'RNA',
  study.name = 'NewStudy',
  experimental.var.name = 'orig.ident',
  reduction = 'umap',
  title_link_text = 'my Lab',
  title_link_url = 'https://github.com/JamesChoi94'
) {
  setup_packages(additional.packages = 'ggplot2')
  if (length(subset.data.paths) != length(subset.data.prefixes)) {
    stop('lengths of subset.data `paths` and `prefixes` do not match')
  }
  full.data <- extract_seurat_data(seurat.object.path = full.data.path)
  
  # Load subset data
  full.data.barcodes <- rownames(full.data$cell.metadata$coldata)
  if (length(subset.data.paths) >= 1) {
    subset.data <- vector(mode = 'list', length = length(subset.data.paths))
    names(subset.data) <- subset.data.prefixes
    for (i in 1:length(subset.data.paths)) {
      subset.data[[i]] <- extract_seurat_data(
        seurat.object.path = subset.data.paths[i],
        include.counts = FALSE,
        default.assay = default.assay
      )
      subset.data.barcodes <- rownames(subset.data[[i]]$cell.metadata$coldata)
      if (!all(subset.data.barcodes %in% full.data.barcodes)) {
        error.message <- paste0('Not all cell barcodes in subset.data `', subset.data.paths[i], '` are contained in the full.data `', full.data.path, '`.')
      }
    }
    merged.coldata <- aggregate_metadata(
      full.data = full.data, 
      subset.data = subset.data
    )
    merged.dimreduc <- aggregate_dimreduc(
      full.data = full.data,
      subset.data = subset.data,
      full.data.prefix = full.data.prefix,
      subset.data.prefixes = subset.data.prefixes
    )
    merged.metadata <- dplyr::full_join(
      x = merged.coldata,
      y = merged.dimreduc,
      by = 'rn'
    )
    rownames(merged.metadata) <- merged.metadata$rn
    merged.metadata$rn <- NULL
  } else {
    merged.metadata <- full.data[['cell.metadata']][['coldata']]
    merged.metadata$rn <- rownames(merged.metadata)
    merged.dimreduc <- full.data[['cell.metadata']][['dimreduc']]
    merged.dimreduc$rn <- rownames(merged.dimreduc)
    merged.metadata <- merge(merged.metadata, merged.dimreduc, by = 'rn')
    rownames(merged.metadata) <- merged.metadata$rn
    merged.metadata$rn <- NULL
  }
  
  dataset.dict <- c(full.data.prefix, subset.data.prefixes)
  names(dataset.dict) <- dataset.dict
  
  is.categorical <- !sapply(merged.metadata, is.numeric)
  cell.groupings <- colnames(merged.metadata)[is.categorical]
  for (i in 1:ncol(merged.metadata)) {
    if (class(merged.metadata[[i]]) == 'character') {
      merged.metadata[[i]] <- factor(merged.metadata[[i]])
    }
  }
  
  all.features <- rownames(full.data$counts)
  
  cell.counts <- as.data.frame(table(merged.coldata[[experimental.var.name]]))
  
  umap_vars <- grep(
    pattern = reduction, 
    x = colnames(merged.metadata),
    value = TRUE,
    ignore.case = TRUE
  )
  label.coords <- merged.metadata[c(cell.groupings, umap_vars)] %>%
    pivot_longer(cols = !contains(match = reduction, ignore.case = TRUE), 
                 names_to = 'variables', 
                 values_to = 'values') %>%
    group_by(variables, values) %>% 
    dplyr::summarise(across(.fns = median)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = 'variables', values_from = 'values') %>% 
    as.data.frame()
  
  window_title = study.name
  title_link_text = title_link_text
  title_link_url = title_link_url
  
  dir.create(path = './data/')
  counts_t <- Matrix::t(full.data$counts)
  HDF5Array::writeTENxMatrix(
    x = counts_t,
    filepath = './data/data.h5',
    verbose = TRUE
  )
  
  rowdata <- full.data$feature.metadata$rowdata
  library.size <- full.data$cell.metadata$library_size
  ncells <- full.data$cell.metadata$ncells
  
  save(merged.metadata,
       rowdata,
       label.coords,
       dataset.dict,
       cell.groupings,
       all.features,
       cell.counts,
       library.size,
       ncells,
       window_title,
       title_link_text,
       title_link_url,
       file = './data/appdata.RData')
}