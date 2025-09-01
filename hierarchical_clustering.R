###############################################################################
# HTT CAG delta — Unsupervised Hierarchical Clustering + ROC report
# Uses first 200 deltas per file; mixture-aware features; Ward.D2 clustering.
###############################################################################

## ---------- User parameters (edit if needed) ----------
input_dir      <- "/PATH/TO/DATA/"
output_pdf     <- "HCA_report.pdf" # change for name of output pdf report
seed           <- 42        
first_n_useful <- 175      # only use the first deltas per dataset
neutral_tau    <- 0.005    # |mean_delta| <= neutral_tau => "neutral"
strong_tau     <- 0.05     # strong pos/neg threshold for features
min_rows       <- 30       # minimum usable deltas after trimming
drop_zeros     <- TRUE     # remove zero-valued deltas?
zero_eps       <- 0        # treat values with |delta| <= zero_eps as zeros

## ------------------------------------------------------

set.seed(seed)

## Packages (auto-install if needed)
install_if_missing <- function(pkgs){
  miss <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
}
install_if_missing(c(
  "readr","dplyr","stringr","purrr","tibble","ggplot2","pROC",
  "dendextend","factoextra","scales","gridExtra","grid","moments",
  "pheatmap","mclust","cluster","tidyr","cowplot"
))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(purrr); library(tibble)
  library(ggplot2); library(pROC); library(dendextend); library(factoextra); library(scales)
  library(gridExtra); library(grid); library(moments); library(pheatmap); library(mclust)
  library(cluster); library(tidyr); library(cowplot)
})

## Helpers
find_delta_col <- function(df){
  nms <- tolower(trimws(names(df)))
  idx <- which(nms == "delta" | stringr::str_detect(nms, "^delta$"))
  if (length(idx)) idx[1] else NA_integer_
}

extract_features <- function(d_in, dataset_name){
  d <- d_in[is.finite(d_in)]
  if (drop_zeros) d <- d[abs(d) > zero_eps]   
  d <- head(d, first_n_useful)
  if (length(d) < min_rows) return(NULL)

  # Robust summaries
  tmn   <- mean(d, trim = 0.1)
  med   <- median(d)
  sdev  <- sd(d)
  madv  <- mad(d, constant = 1.4826, na.rm = TRUE)
  skew  <- if (length(unique(d)) >= 3) moments::skewness(d) else 0
  iqr_  <- IQR(d)
  ppos  <- mean(d > 0)
  pneg  <- mean(d < 0)
  pstrp <- mean(d >=  strong_tau)
  pstrn <- mean(d <= -strong_tau)
  pneut <- mean(abs(d) <= neutral_tau)
  mean_abs <- mean(abs(d))

  # Mixture features (1–2 comps), safe on constants
  mix_k <- 1L; mix_mu1 <- mean(d); mix_mu2 <- NA_real_
  mix_pi1 <- 1; mix_pi2 <- NA_real_
  mix_sep <- 0; mix_entropy <- 0
  if (length(unique(d)) >= 2 && is.finite(sdev) && sdev > 0) {
    fit <- tryCatch(Mclust(d, G = 1:2, modelNames = c("E","V"), verbose = FALSE), error = function(e) NULL)
    if (!is.null(fit)) {
      mix_k <- fit$G
      ord   <- order(fit$parameters$mean)
      mu    <- as.numeric(fit$parameters$mean)[ord]
      piw   <- as.numeric(fit$parameters$pro)[ord]
      if (mix_k == 1) {
        mix_mu1 <- mu[1]; mix_pi1 <- 1
      } else {
        mix_mu1 <- mu[1]; mix_mu2 <- mu[2]
        mix_pi1 <- piw[1]; mix_pi2 <- piw[2]
        pooled_sd <- sd(d)
        mix_sep <- if (is.finite(pooled_sd) && pooled_sd > 0) abs(mu[2]-mu[1]) / pooled_sd else 0
        w <- pmax(piw, 1e-12); mix_entropy <- -sum(w * log(w))
      }
    }
  }

  tibble(
    dataset     = dataset_name,
    n           = length(d),
    mean_delta  = mean(d),
    trimmed_mean= tmn,
    median_delta= med,
    sd_delta    = sdev,
    mad_delta   = madv,
    skew_delta  = skew,
    iqr_delta   = iqr_,
    mean_abs    = mean_abs,
    prop_pos    = ppos,
    prop_neg    = pneg,
    prop_neutral= pneut,
    prop_strong_pos = pstrp,
    prop_strong_neg = pstrn,
    mix_k       = mix_k,
    mix_mu1     = mix_mu1,
    mix_mu2     = mix_mu2,
    mix_pi1     = mix_pi1,
    mix_pi2     = mix_pi2,
    mix_sep     = mix_sep,
    mix_entropy = mix_entropy
  )
}

## 1) Ingest CSVs, keep first N deltas, compute features
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No CSV files found in: ", input_dir)

raw_for_plots <- list()
features_list <- list()
for (fp in files) {
  nm <- tools::file_path_sans_ext(basename(fp))
  res <- tryCatch({
    dat <- suppressMessages(readr::read_csv(fp, show_col_types = FALSE))
    ci  <- find_delta_col(dat); if (is.na(ci)) stop("no 'delta' column")
del <- suppressWarnings(as.numeric(dat[[ci]]))
del <- del[is.finite(del)]
if (drop_zeros) del <- del[abs(del) > zero_eps]   # <-- NEW

if (!length(del)) stop("no numeric deltas after filtering")

trimmed <- head(del, first_n_useful)              # now: first 200 *non-zero* deltas
if (length(trimmed) < min_rows) stop("too few deltas after trimming")

    raw_for_plots[[nm]] <- trimmed
    extract_features(del, nm)
  }, error = function(e){ warning(sprintf("Skipping %s: %s", basename(fp), conditionMessage(e))); NULL })
  if (!is.null(res)) features_list[[nm]] <- res
}
if (!length(features_list)) stop("No usable datasets with 'delta' values after trimming.")
features <- bind_rows(features_list) %>% arrange(dataset)

## 2) Coarse “truth” (for evaluation only; training stays unsupervised)
features <- features %>%
  mutate(true_label = case_when(
    mean_delta >  neutral_tau ~ "expansion",
    mean_delta < -neutral_tau ~ "contraction",
    TRUE                      ~ "neutral"
  ))

## 3) Train/Test split (dataset-level 80/20)
ids <- features$dataset; n_total <- length(ids); n_train <- max(1, floor(0.8*n_total))
ids_exp <- features %>% filter(true_label=="expansion") %>% pull(dataset)
ids_con <- features %>% filter(true_label=="contraction") %>% pull(dataset)
ids_neu <- features %>% filter(true_label=="neutral")    %>% pull(dataset)
take <- function(v,k){ if (!length(v)) character(0) else sample(v, min(length(v), k)) }
prop_exp <- length(ids_exp)/n_total; prop_con <- length(ids_con)/n_total
tr_exp <- take(ids_exp, round(n_train*prop_exp)); tr_con <- take(ids_con, round(n_train*prop_con))
tr_neu <- take(ids_neu, n_train - length(tr_exp) - length(tr_con))
train_ids <- unique(c(tr_exp, tr_con, tr_neu))
if (length(train_ids) < n_train) train_ids <- c(train_ids, sample(setdiff(ids, train_ids), n_train - length(train_ids)))
test_ids  <- setdiff(ids, train_ids)
train_df <- features %>% filter(dataset %in% train_ids)
test_df  <- features %>% filter(dataset %in% test_ids)

## 4) Build feature matrices, impute, scale (train stats only)
feature_cols <- c(
  "mean_delta","trimmed_mean","median_delta","sd_delta","mad_delta","skew_delta",
  "iqr_delta","mean_abs","prop_pos","prop_neg","prop_neutral",
  "prop_strong_pos","prop_strong_neg","mix_k","mix_mu1","mix_mu2",
  "mix_pi1","mix_pi2","mix_sep","mix_entropy"
)
feature_cols <- intersect(feature_cols, names(features))
if (length(feature_cols) < 2) stop("Too few feature columns present.")

X_all   <- as.matrix(features[,  feature_cols, drop = FALSE])
X_train <- as.matrix(train_df[,  feature_cols, drop = FALSE])
X_test  <- as.matrix(test_df[,   feature_cols, drop = FALSE])

impute_median <- function(M){
  for (j in seq_len(ncol(M))) {
    cj <- M[,j]; med <- suppressWarnings(stats::median(cj, na.rm = TRUE)); if (!is.finite(med)) med <- 0
    cj[is.na(cj) | !is.finite(cj)] <- med; M[,j] <- cj
  }; M
}
X_all   <- impute_median(X_all)
X_train <- impute_median(X_train)
X_test  <- impute_median(X_test)

rownames(X_all)   <- features$dataset
rownames(X_train) <- train_df$dataset
rownames(X_test)  <- test_df$dataset

sc_center <- colMeans(X_train, na.rm = TRUE)
sc_scale  <- apply(X_train, 2, sd); sc_scale[!is.finite(sc_scale) | sc_scale == 0] <- 1
scale_mat <- function(M) sweep(sweep(M, 2, sc_center, "-"), 2, sc_scale, "/")
Z_train <- scale_mat(X_train); Z_test <- scale_mat(X_test); Z_all <- scale_mat(X_all)

## 5) HCA on training; auto-label clusters
dist_train <- dist(Z_train, method = "euclidean")
hc_train   <- hclust(dist_train, method = "ward.D2")
cl_train   <- cutree(hc_train, k = 2)

mean_by_cl <- train_df %>% mutate(cluster = cl_train) %>%
  group_by(cluster) %>% summarise(avg_mean_delta = mean(mean_delta), .groups="drop") %>% arrange(desc(avg_mean_delta))
cl_names <- setNames(ifelse(rank(-mean_by_cl$avg_mean_delta)==1,"expansion","contraction"), mean_by_cl$cluster)
train_cluster_name <- cl_names[as.character(cl_train)]

# Centroids in standardized space
centroids <- lapply(names(cl_names), function(ci){
  rows <- which(cl_train == as.integer(ci)); colMeans(Z_train[rows,,drop=FALSE])
})
names(centroids) <- names(cl_names); C <- do.call(rbind, centroids)
row_to_name <- cl_names[rownames(C)]
ci_exp <- which(row_to_name == "expansion"); ci_con <- which(row_to_name == "contraction")

## 6) Assign test by nearest centroid; compute ROC (non-neutral only)
dist_to_C <- function(z, Cmat) apply(Cmat, 1, function(crow) sqrt(sum((z - crow)^2)))
assign_test <- t(apply(Z_test, 1, dist_to_C, Cmat = C)) %>% as.data.frame()
score <- assign_test[[ci_con]] - assign_test[[ci_exp]]    # higher => more expansion-like
pred_label_test <- ifelse(assign_test[[ci_exp]] < assign_test[[ci_con]], "expansion","contraction")
test_eval <- test_df %>% mutate(pred_cluster = pred_label_test, score = score)

roc_df <- test_eval %>% filter(true_label %in% c("expansion","contraction")) %>%
  mutate(y = ifelse(true_label=="expansion",1,0))
roc_obj <- NULL; auc_val <- NA
if (nrow(roc_df) >= 2 && length(unique(roc_df$y))==2) {
  roc_obj <- pROC::roc(response = roc_df$y, predictor = roc_df$score, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
}

## 7) PDF report
pdf(output_pdf, width = 9, height = 7)

# Title
grid.newpage()
title <- "HCA_report.pdf - HCA of HTT CAG delta datasets"
subt  <- paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
                "   |   Folder: ", normalizePath(input_dir, mustWork = FALSE))
txt   <- paste0(
  "Datasets: ", nrow(features), " (first ", first_n_useful, " deltas each)\n",
  "Train/Test: ", length(train_ids), "/", length(test_ids), "  |  Neutral band: ±", neutral_tau, "\n",
  "Clustering: Ward.D2 on standardized robust + mixture features\n"
)
grid.draw(textGrob(title, gp=gpar(fontsize=18, fontface="bold"), y=0.9))
grid.draw(textGrob(subt,  gp=gpar(fontsize=10, col="grey30"), y=0.85))
grid.draw(textGrob(txt,   gp=gpar(fontsize=12), y=0.6))

# Dendrogram (training), colored labels
dend <- as.dendrogram(hc_train)
labels(dend) <- rownames(Z_train)
dend <- dendextend::color_branches(dend, k = 2) # dendrogram threshold
lab_cols <- setNames(ifelse(train_cluster_name == "expansion","steelblue","firebrick"), rownames(Z_train))
dend <- dendextend::set(dend, "labels_col", value = lab_cols[labels(dend)])
plot(dend, main = "Training dendrogram (Ward.D2)\nLabels: expansion (blue) vs contraction (red)")

# Silhouette (if feasible)
if (nrow(Z_train) >= 3 && length(unique(cl_train)) >= 2) {
  sil <- cluster::silhouette(cl_train, dist_train)
  print(factoextra::fviz_silhouette(sil) + ggtitle("Silhouette (training)"))
} else {
  plot.new(); title(main = "Silhouette (training)")
  grid.text("Not enough training datasets for silhouette.", y = 0.5)
}

# Heatmap of standardized features (all datasets)
if (nrow(Z_all) >= 2) {
  dist_all <- dist(Z_all, method = "euclidean")
  hc_all   <- hclust(dist_all, method = "ward.D2")
tmp <- features
tmp$split <- ifelse(tmp$dataset %in% train_ids, "train", "test")
tmp$label <- paste0(tmp$split, " / ", tmp$true_label)
row_annot <- data.frame(label = tmp$label, row.names = tmp$dataset)

  ann_df <- data.frame(SplitLabel = row_annot[rownames(Z_all), "label", drop = TRUE])
  rownames(ann_df) <- rownames(Z_all)
  pheatmap::pheatmap(
    Z_all, cluster_rows = hc_all, cluster_cols = TRUE,
    main = "Standardized feature heatmap (all datasets)",
    annotation_row = ann_df, show_rownames = TRUE, fontsize_row = 7
  )
} else {
  plot.new(); title(main = "Standardized feature heatmap")
  grid.text("Need at least 2 datasets.", y = 0.5)
}

# ROC curve (test set, non-neutral)
if (!is.null(roc_obj)) {
  plot(roc_obj, main = sprintf("ROC (test, excluding neutral) — AUC = %.3f", auc_val))
} else {
  plot.new(); title(main = "ROC (test set)")
  grid.text("Not enough non-neutral test datasets to compute ROC.", y = 0.5)
}

# Distribution of mean deltas
p1 <- ggplot(features %>% mutate(split = ifelse(dataset %in% train_ids,"train","test")),
             aes(x = mean_delta, fill = split)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20) +
  geom_vline(xintercept = c(-neutral_tau, neutral_tau), linetype = "dashed") +
  labs(title = "Distribution of dataset mean deltas (first N only)",
       subtitle = paste0("Neutral band: ±", neutral_tau), x = "mean(delta)", y = "Count") +
  theme_minimal()
print(p1)

# Small-multiples: per-dataset delta histograms (up to 12)
if (length(raw_for_plots)) {
  show_names <- head(sort(names(raw_for_plots)), 12)
  plots <- lapply(show_names, function(nm){
    df <- tibble(val = raw_for_plots[[nm]])
    ggplot(df, aes(x = val)) +
      geom_histogram(bins = 40, alpha = 0.7) +
      geom_vline(xintercept = c(-neutral_tau, neutral_tau), linetype = "dashed") +
      labs(title = nm, x = "delta (first N)", y = "count") +
      theme_minimal(base_size = 9)
  })
  grid.newpage(); title(main = "Per-dataset delta histograms (first N)")
  print(cowplot::plot_grid(plotlist = plots, ncol = 3))
}

# Test set summary (non-neutral only)
conf_like <- test_eval %>%
  filter(true_label %in% c("expansion","contraction")) %>%
  count(true_label, pred_cluster, name = "n") %>%
  tidyr::complete(true_label = c("contraction","expansion"),
                  pred_cluster = c("contraction","expansion"),
                  fill = list(n = 0)) %>%
  arrange(true_label, pred_cluster)
plot.new(); title(main = "Test set summary (non-neutral only)")
gridExtra::grid.table(conf_like)

# Session info
plot.new(); title(main = "Session info")
si <- capture.output(sessionInfo())
grid.text(paste(si, collapse = "\n"), x=0.02, y=0.98, just=c("left","top"), gp=gpar(fontsize=7))

dev.off()

## 8) Save features/predictions CSV + console pointers
# Add split without dplyr
features_split <- within(features, split <- ifelse(dataset %in% train_ids, "train", "test"))

# Subset test_eval safely (only keep cols that exist)
keep_cols <- intersect(c("dataset", "pred_cluster", "score"), names(test_eval))
test_keep <- if (length(keep_cols)) test_eval[, keep_cols, drop = FALSE] else
             data.frame(dataset = character(), pred_cluster = character(), score = numeric())

# Left join using base merge (all.x=TRUE keeps all features rows)
out_tbl <- merge(features_split, test_keep, by = "dataset", all.x = TRUE, sort = FALSE)

readr::write_csv(out_tbl, "HCA_per_dataset_features_and_results.csv")


cat("\nReport: ", normalizePath(output_pdf, mustWork = FALSE), "\n", sep = "")
cat("Per-dataset table: HCA_per_dataset_features_and_results.csv\n")
cat("Train/Test sizes: ", length(train_ids), "/", length(test_ids), "\n", sep = "")
if (!is.null(roc_obj)) cat(sprintf("Test ROC AUC (non-neutral): %.3f\n", auc_val))
