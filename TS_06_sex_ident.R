library(dplyr)
library(readr)

# Read pseudobulk files

read_all_rds_to_env <- function(path, envir = .GlobalEnv) {
  files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
  for (f in files) {
    obj_name <- tools::file_path_sans_ext(basename(f))
    assign(obj_name, readRDS(f), envir = envir)
  }
}

read_all_rds_to_env(path = "Data/pseudobulk_per_major_cell_type")

projids_used <- as.character(unique(readRDS("projids_used.rds")))
meta <- read.csv("Data/Metadata.csv")
meta$projid <- as.character(meta$projid)

meta <- dplyr::semi_join(meta, data.frame(projid = projids_used), by = "projid")
meta <- meta[meta$projid %in% projids_used, ]


plot_xist_vs_sry <- function(cell_type,
                             envir = .GlobalEnv,
                             gene_x = "XIST",
                             gene_y = "SRY",
                             transform = c("log1p","none"),
                             point_alpha = 0.8,
                             label_points = FALSE,
                             label_mixed_only = TRUE,
                             mixed_radius = 0.5,
                             min_opposite_neighbors = 5,
                             label_isolated = FALSE,
                             iso_radius = 0.35,
                             alias_map = list(XIST = c("Xist"),
                                              SRY  = c("TDF")),
                             meta_df = NULL,
                             meta_name = "meta",
                             sex_col = "msex",
                             sex_map = NULL,
                             sex_palette = c(Female = "#ff8fab", Male = "#219ebc")) {
  
  transform <- match.arg(transform)
  
  obj_name <- if (exists(cell_type, envir = envir, inherits = FALSE)) {
    cell_type
  } else if (exists(paste0(cell_type, "_pseudo"), envir = envir, inherits = FALSE)) {
    paste0(cell_type, "_pseudo")
  } else {
    stop("Could not find object '", cell_type,
         "' or '", paste0(cell_type, "_pseudo"), "' in the given environment.", call. = FALSE)
  }
  
  mat <- get(obj_name, envir = envir)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have rownames (genes) and colnames (projid).", call. = FALSE)
  }
  
  resolve_gene <- function(g, aliases = NULL) {
    rn <- rownames(mat)
    hit <- which(toupper(rn) == toupper(g))
    if (length(hit) >= 1) return(list(idx = hit[1], used = rn[hit[1]]))
    if (!is.null(aliases) && length(aliases)) {
      for (a in aliases) {
        hit <- which(toupper(rn) == toupper(a))
        if (length(hit) >= 1) return(list(idx = hit[1], used = rn[hit[1]]))
      }
    }
    sugg <- unique(c(
      grep(g, rn, ignore.case = TRUE, value = TRUE),
      tryCatch(agrep(g, rn, ignore.case = TRUE, value = TRUE, max.distance = 0.2),
               error = function(e) character(0))
    ))
    if (length(sugg) > 8) sugg <- sugg[1:8]
    looks_ensg   <- any(grepl("^ENSG\\d+", rn))
    looks_symbol <- any(grepl("^[A-Za-z0-9]+$", rn))
    hint <- if (looks_ensg && !looks_symbol) {
      "Row names look like Ensembl IDs (e.g., ENSG...). If you supplied a symbol, map symbols → Ensembl first."
    } else if (!looks_ensg && looks_symbol) {
      "Row names look like gene symbols. Check species/case/alias."
    } else {
      "Check whether row names are symbols or Ensembl IDs."
    }
    stop(sprintf("Gene '%s' not found in %s.\n%s\nClosest matches: %s",
                 g, obj_name, hint,
                 if (length(sugg)) paste(sugg, collapse = ", ") else "none"),
         call. = FALSE)
  }
  
  gx <- resolve_gene(gene_x, aliases = alias_map[[gene_x]])
  gy <- resolve_gene(gene_y, aliases = alias_map[[gene_y]])
  
  x <- as.numeric(mat[gx$idx, ])
  y <- as.numeric(mat[gy$idx, ])
  projid <- colnames(mat)
  
  if (transform == "log1p") {
    x <- log1p(x); y <- log1p(y)
    xlab <- paste0(gx$used, " (log1p)")
    ylab <- paste0(gy$used, " (log1p)")
  } else {
    xlab <- gx$used; ylab <- gy$used
  }
  
  df <- data.frame(projid = as.character(projid), x = x, y = y, stringsAsFactors = FALSE)
  
  if (is.null(meta_df)) {
    if (!exists(meta_name, envir = envir, inherits = FALSE)) {
      stop("meta_df not provided and object '", meta_name, "' not found in 'envir'.", call. = FALSE)
    }
    meta_df <- get(meta_name, envir = envir)
  }
  if (!all(c("projid", sex_col) %in% colnames(meta_df))) {
    stop("meta_df must contain columns 'projid' and '", sex_col, "'.", call. = FALSE)
  }
  
  meta_sub <- meta_df[, c("projid", sex_col)]
  meta_sub$projid <- as.character(meta_sub$projid)
  
  if (is.null(sex_map)) {
    raw_vals <- unique(na.omit(meta_sub[[sex_col]]))
    if (all(raw_vals %in% c(0, 1))) {
      sex_map <- c("0" = "Female", "1" = "Male")
    } else if (all(raw_vals %in% c(1, 2))) {
      sex_map <- c("1" = "Male", "2" = "Female")
    } else {
      sex_map <- c("F" = "Female","f" = "Female","Female" = "Female",
                   "M" = "Male",  "m" = "Male",  "Male"   = "Male")
    }
  }
  
  meta_sub$Sex <- unname(sex_map[as.character(meta_sub[[sex_col]])])
  meta_sub$Sex[is.na(meta_sub$Sex)] <- NA_character_
  
  df <- dplyr::inner_join(df, meta_sub[, c("projid", "Sex")], by = "projid")
  df$Sex <- factor(df$Sex, levels = c("Female", "Male"))
  
  # ----------------- LABELING RULES -----------------
  n <- nrow(df)
  df$to_label <- FALSE
  
  if (label_points && label_mixed_only) {
    idxF <- which(df$Sex == "Female")
    idxM <- which(df$Sex == "Male")
    if (length(idxF) > 0 && length(idxM) > 0) {
      dFM <- sqrt(outer(df$x[idxF], df$x[idxM], "-")^2 +
                    outer(df$y[idxF], df$y[idxM], "-")^2)
      cntF <- if (length(idxM) > 0) rowSums(dFM <= mixed_radius, na.rm = TRUE) else rep(0L, length(idxF))
      cntM <- if (length(idxF) > 0) colSums(dFM <= mixed_radius, na.rm = TRUE) else rep(0L, length(idxM))
      df$to_label[idxF] <- cntF >= min_opposite_neighbors
      df$to_label[idxM] <- cntM >= min_opposite_neighbors
    }
  } else if (label_points && !label_mixed_only) {
    df$to_label <- TRUE
  }
  
  if (label_isolated) {
    if (n <= 1) {
      iso_flag <- rep(TRUE, n)
    } else {
      dx <- outer(df$x, df$x, "-")
      dy <- outer(df$y, df$y, "-")
      dist_mat <- sqrt(dx^2 + dy^2)
      diag(dist_mat) <- Inf  
      nb_cnt <- rowSums(dist_mat <= iso_radius, na.rm = TRUE)
      iso_flag <- nb_cnt == 0L
    }
    df$to_label <- df$to_label | iso_flag
  }
  
  # ----------------- PLOT -----------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.", call. = FALSE)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = Sex)) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::scale_color_manual(values = sex_palette, drop = FALSE, na.translate = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0(obj_name, ": ", gx$used, " vs ", gy$used, " per projid"),
      x = xlab, y = ylab, color = "Sex"
    )
  
  # labels
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Please install 'ggrepel' for labels.", call. = FALSE)
  if (any(df$to_label)) {
    p <- p + ggrepel::geom_text_repel(
      data = subset(df, to_label),
      ggplot2::aes(label = projid),
      size = 5, max.overlaps = 100,
      show.legend = FALSE
    )
  }
  
  return(p)
}


# usage
p <- plot_xist_vs_sry("Exc",
                      transform = "log1p",
                      label_points = TRUE,
                      label_mixed_only = TRUE,
                      mixed_radius = 0.1,
                      min_opposite_neighbors = 2,
                      label_isolated = TRUE,
                      iso_radius = 0.1,
                      gene_y = "UTY")

library(ggplot2)

p <- p + theme(text = element_text(family = "Arial"))

ggsave("plot_xist_vs_uty.pdf", plot = p,
       width = 6, height = 5, units = "in",
       device = cairo_pdf)
ggsave("plot_xist_vs_uty.svg", plot = p,
       width = 6, height = 5, units = "in")



plot_xist_vs_sry("Exc",
                 transform = "log1p",
                 label_points = TRUE,
                 label_mixed_only = TRUE,
                 mixed_radius = 0.4,
                 min_opposite_neighbors = 5,
                 label_isolated = TRUE,   
                 iso_radius = 0.1,      
                 gene_y = "DDX3Y")

plot_xist_vs_sry("Inh",
                 transform = "log1p",
                 label_points = TRUE,
                 label_mixed_only = TRUE,
                 mixed_radius = 0.1,
                 min_opposite_neighbors = 5,
                 label_isolated = TRUE,  
                 iso_radius = 15,     
                 gene_y = "DDX3Y",
                 gene_x = "KDM6A")