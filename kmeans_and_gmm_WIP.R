# ========================= USER INPUT =========================
train_folder <- "/Users/angusdixon/masters/met591_diss/all_data/ML0708/training_data"   # folder of CSVs (first col = repeat length, col 'delta') (currently 8/11 datasets)
test_folder  <- "/Users/angusdixon/masters/met591_diss/all_data/ML0708/test_data"       # folder of CSVs for evaluation; can equal train_folder (currently 3/11 datasets)

method        <- "gmm"            # "gmm" (mclust) or "kmeans"
K             <- 3                # used if auto_select_K = FALSE
auto_select_K <- TRUE             # for GMM: let mclust choose K by BIC (range below) (k-means = TRUE)

bin_width <- 1                    # repeat-length bin size (k-means = 1)
tau       <- 0.0075               # proxy threshold for Δ labels (bin-level ROC) (k-means = 0.001)
smooth_k  <- 5                    # small pre-smoothing for plots only (rolling mean) (k-means = 5)

set.seed(42)

# ========================= SETUP ==============================
pkgs <- c("tidyverse","ggplot2","zoo","pROC","patchwork","scales","readr","mclust")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE)) # install packages

stopifnot(dir.exists(train_folder)) # stop if test or train folder doesnt exist
stopifnot(dir.exists(test_folder))

train_files <- list.files(train_folder, pattern="\\.csv$", full.names=TRUE)
test_files  <- list.files(test_folder,  pattern="\\.csv$", full.names=TRUE)
stopifnot(length(train_files) > 0, length(test_files) > 0)

out_dir <- file.path(getwd(), "unsup_model_out")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ========================= HELPERS ============================
# confirm presence and format delta column
read_delta <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  if (!"delta" %in% names(df)) stop("File lacks a 'delta' column: ", path)
  xname <- names(df)[1]
  df <- df %>% arrange(.data[[xname]])
  tibble(x = as.numeric(df[[xname]]), delta = as.numeric(df$delta))
}

roll_smooth <- function(y, k=3) {
  ys <- zoo::rollmean(y, k = k, align = "center", fill = NA)
  ifelse(is.na(ys), y, ys)
}

# Build common bins from TRAIN data, then apply to all
compute_bins <- function(files, bin_width) {
  xr <- purrr::map(files, ~read_delta(.x)$x) %>% unlist()
  xmin <- floor(min(xr, na.rm=TRUE))
  xmax <- ceiling(max(xr, na.rm=TRUE))
  if (!is.finite(xmin) || !is.finite(xmax) || xmax <= xmin) stop("Invalid x range in training data.")
  breaks <- seq(xmin, xmax + bin_width, by = bin_width)
  list(breaks = breaks, mids = head(breaks, -1) + diff(breaks)/2)
}

# Bin a single series into mean Δ per bin
bin_series <- function(x, y, breaks) {
  cut_idx <- cut(x, breaks = breaks, include.lowest = TRUE, right = FALSE)
  as.numeric(tapply(y, cut_idx, function(v) mean(v, na.rm=TRUE)))
}

# Feature for one file (vector of bin-mean Δ, NA -> 0)
file_to_features <- function(path, breaks) {
  dat <- read_delta(path)
  fv <- bin_series(dat$x, dat$delta, breaks)
  fv[is.na(fv)] <- 0
  fv
}

# ROC plotting (per-sample)
roc_plot_from <- function(labels01, scores, title_txt) {
  ok <- !is.na(labels01) & !is.na(scores)
  lbl <- labels01[ok]; sc <- scores[ok]

  if (length(lbl) < 2 || length(unique(lbl)) < 2) {
    return(
      ggplot() +
        annotate("text", x=0.5, y=0.5, label="ROC unavailable (one class)", size=4) +
        xlim(0,1) + ylim(0,1) +
        labs(title = title_txt, x="False Positive Rate", y="True Positive Rate") +
        theme_minimal()
    )
  }
  roc_obj <- pROC::roc(lbl, sc, quiet=TRUE, levels=c(0,1), direction=">")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_auc  <- tryCatch(pROC::ci.auc(roc_obj), error=function(e) NA)
  cr <- pROC::coords(roc_obj, x="all", ret=c("specificity","sensitivity"), transpose=FALSE) |> as.data.frame()
  df_roc <- tibble(fpr = 1 - as.numeric(cr$specificity),
                   tpr =       as.numeric(cr$sensitivity)) %>% dplyr::filter(stats::complete.cases(fpr,tpr))
  if (nrow(df_roc) < 2) {
    return(
      ggplot() +
        annotate("text", x=0.5, y=0.5, label=sprintf("AUC = %.3f (degenerate)", auc_val), size=4) +
        xlim(0,1) + ylim(0,1) +
        labs(title = title_txt, x="False Positive Rate", y="True Positive Rate") +
        theme_minimal()
    )
  }
  ann_txt <- if (all(!is.na(ci_auc))) sprintf("AUC = %.3f (95%% CI %.3f–%.3f)", auc_val, ci_auc[1], ci_auc[3])
             else sprintf("AUC = %.3f", auc_val)

  ggplot(df_roc, aes(fpr, tpr)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    geom_path() +
    annotate("text", x=0.65, y=0.1, hjust=0, label=ann_txt) +
    labs(title=title_txt, x="False Positive Rate", y="True Positive Rate") +
    coord_equal() + theme_minimal()
}

# Minimal overall-ROC plot (single page)
overall_roc_min_plot <- function(labels01, scores, title_txt) {
  ok <- !is.na(labels01) & !is.na(scores)
  lbl <- labels01[ok]; sc <- scores[ok]

  if (length(lbl) < 2 || length(unique(lbl)) < 2) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "Overall ROC unavailable (one class)", size = 4, color = "black") +
        xlim(0,1) + ylim(0,1) +
        labs(title = title_txt, x = "False Positive Rate", y = "True Positive Rate") +
        theme_classic(base_size = 12) +
        theme(panel.grid = element_blank(),
              axis.text = element_text(color = "black"),
              axis.title = element_text(color = "black"),
              plot.title = element_text(color = "black"))
    )
  }

  roc_obj <- pROC::roc(lbl, sc, quiet = TRUE, levels = c(0,1), direction = ">")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_auc  <- tryCatch(pROC::ci.auc(roc_obj), error = function(e) NA)

  cr <- pROC::coords(roc_obj, x = "all",
                     ret = c("specificity","sensitivity"), transpose = FALSE) |>
        as.data.frame()

  df_roc <- tibble::tibble(
    fpr = 1 - as.numeric(cr$specificity),
    tpr =       as.numeric(cr$sensitivity)
  ) |> dplyr::filter(stats::complete.cases(fpr, tpr))

  auc_txt <- if (all(!is.na(ci_auc))) {
    sprintf("AUC = %.3f\n95%% CI %.3f–%.3f", auc_val, ci_auc[1], ci_auc[3])
  } else sprintf("AUC = %.3f", auc_val)

  ggplot(df_roc, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    geom_path(linewidth = 1, color = "red") +                # red ROC line
    annotate("label", x = 0.65, y = 0.2, label = auc_txt,    # simple corner box
             color = "black", fill = "white", size = 3, hjust = 0, vjust = 0) +
    labs(title = title_txt, x = "False Positive Rate", y = "True Positive Rate") +
    coord_equal() +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          plot.title = element_text(color = "black"))
}

# Human-readable cluster label from center mean (used to tag cluster types)
label_cluster_type <- function(center_vec, tau) {
  m <- mean(center_vec, na.rm=TRUE)
  if (m >  +tau) "expansion-dominant"
  else if (m < -tau) "contraction-dominant"
  else "neutral-like"
}

# ========================= TRAINING DATA -> FEATURES =======================
bins   <- compute_bins(train_files, bin_width = bin_width)
breaks <- bins$breaks
mids   <- bins$mids
p      <- length(mids)
if (p < 2) stop("Not enough bins/features (p < 2). Increase `bin_width` or ensure x-range spans multiple bins.")

X_train <- purrr::map(train_files, ~file_to_features(.x, breaks)) |> do.call(what = rbind)
colnames(X_train) <- paste0("bin_", seq_len(p))
rownames(X_train) <- basename(train_files)

# Coerce to numeric matrix
X_train <- as.matrix(X_train); storage.mode(X_train) <- "double"

# Standardize features (per-bin z-score) using TRAIN only
train_mean <- colMeans(X_train, na.rm = TRUE)
train_sd   <- apply(X_train, 2, sd, na.rm = TRUE); train_sd[train_sd == 0] <- 1

X_train_s  <- scale(X_train, center = train_mean, scale = train_sd)
X_train_s  <- as.matrix(X_train_s); storage.mode(X_train_s) <- "double"

n_train <- nrow(X_train_s)

# ========================= FIT MODEL (K-means or GMM) ======================
model <- list(method = method)

if (tolower(method) == "kmeans") {
  K <- max(2, min(K, n_train))  # safe range
  km <- kmeans(X_train_s, centers = K, nstart = 50, iter.max = 1000)
  centers_s    <- km$centers                       # K x p
  centers_orig <- sweep(centers_s, 2, train_sd, `*`)
  centers_orig <- sweep(centers_orig, 2, train_mean, `+`)
  cluster_types <- apply(centers_orig, 1, label_cluster_type, tau = tau)

  model$kmeans        <- km
  model$centers_s     <- centers_s
  model$centers_orig  <- centers_orig
  model$cluster_types <- cluster_types

} else {
  # ----- GMM via mclust -----
  if (auto_select_K) {
    maxK <- max(2, min(9, n_train - 1))
    mc <- Mclust(data = X_train_s, G = 1:maxK)
    K  <- mc$G
  } else {
    K <- max(2, min(K, n_train - 1))
    mc <- Mclust(data = X_train_s, G = K)
  }

  # Means shape handling
  mu <- mc$parameters$mean
  if (is.null(dim(mu))) {
    centers_s <- matrix(mu, nrow = K, ncol = 1, byrow = FALSE) # p=1 case
  } else {
    centers_s <- t(mu)  # K x p (mclust stores means as p x K)
  }

  centers_orig <- sweep(centers_s, 2, train_sd, `*`)
  centers_orig <- sweep(centers_orig, 2, train_mean, `+`)
  cluster_types <- apply(centers_orig, 1, label_cluster_type, tau = tau)

  model$mclust        <- mc
  model$centers_s     <- centers_s
  model$centers_orig  <- centers_orig
  model$cluster_types <- cluster_types
}

model$breaks     <- breaks
model$mids       <- mids
model$train_mean <- train_mean
model$train_sd   <- train_sd
model$bin_width  <- bin_width
model$tau        <- tau
model$K          <- K

model_path <- file.path(out_dir, paste0("unsup_model_", method, "_K", K, ".rds"))
saveRDS(model, model_path)
message("Saved model to: ", model_path)

# ========================= TRAIN SUMMARY CSVs ==============================
if (tolower(method) == "kmeans") {
  train_cls <- model$kmeans$cluster
  # soft membership via inverse-distance weights (for reporting only)
  train_post <- t(sapply(1:nrow(X_train_s), function(i) {
    d <- apply(model$centers_s, 1, function(cn) sqrt(sum((X_train_s[i,] - cn)^2)))
    w <- 1/(d^2 + 1e-6); w/sum(w)
  }))
} else {
  pred_train <- predict(model$mclust, newdata = X_train_s)
  if (is.null(pred_train$z)) stop("GMM training failed: null posteriors. Try different K or check inputs.")
  train_cls  <- pred_train$classification
  train_post <- pred_train$z
}

train_sum <- tibble(file = rownames(X_train_s),
                    cluster = as.integer(train_cls),
                    cluster_type = model$cluster_types[train_cls]) %>%
             bind_cols(as_tibble(train_post, .name_repair = ~paste0("pC", seq_len(ncol(train_post)))))

write_csv(train_sum, file.path(out_dir, "train_assignments.csv"))
write_csv(as_tibble(model$centers_orig) %>% mutate(cluster = row_number(),
           cluster_type = model$cluster_types),
          file.path(out_dir, "cluster_centers.csv"))

# ========================= TEST EVALUATION & REPORT =======================
pdf_fn <- file.path(out_dir, "test_report.pdf")
grDevices::pdf(pdf_fn, width = 11, height = 8.5)
on.exit(grDevices::dev.off(), add = TRUE)

# Collect labels/scores across ALL test files for a global ROC
all_lbl_exp <- integer(0); all_scr_exp <- numeric(0)
all_lbl_con <- integer(0); all_scr_con <- numeric(0)

metrics <- tibble()
pred_bins_all <- tibble()

predict_from_model <- function(fv_s, model) {
  # fv_s: numeric vector (scaled), length p
  if (tolower(model$method) == "kmeans") {
    dists <- apply(model$centers_s, 1, function(cn) sqrt(sum((fv_s - cn)^2)))
    cls <- which.min(dists)
    pred <- model$centers_orig[cls, ]
    probs <- {w <- 1/(dists^2 + 1e-6); w/sum(w)}
    list(pred_bin_delta = as.numeric(pred), cluster = cls, cluster_probs = probs)
  } else {
    fv_s_mat <- matrix(fv_s, nrow = 1)
    pr <- predict(model$mclust, newdata = fv_s_mat)
    if (is.null(pr$z)) stop("GMM predict failed: null posteriors for a test sample.")
    post <- drop(pr$z)              # length K
    cls  <- pr$classification[1]
    pred <- as.numeric(drop(post %*% model$centers_orig))
    list(pred_bin_delta = pred, cluster = cls, cluster_probs = post)
  }
}

for (f in test_files) {
  bn  <- basename(f)
  dat <- read_delta(f)

  # Features -> scale with TRAIN stats
  fv   <- file_to_features(f, breaks = model$breaks)
  fv_s <- as.numeric((fv - model$train_mean) / model$train_sd)
  if (any(!is.finite(fv_s))) fv_s[!is.finite(fv_s)] <- 0
  stopifnot(length(fv_s) == length(model$mids))

  # Predict expected Δ per bin from model
  pred   <- predict_from_model(fv_s, model)
  pred_bins <- pred$pred_bin_delta
  cls       <- as.integer(pred$cluster)
  probs     <- as.numeric(pred$cluster_probs)
  cls_type  <- model$cluster_types[cls]

  # Proxy labels at BIN level for ROC
  y_bin <- fv
  y_lbl_exp <- ifelse(y_bin >  model$tau, 1L, ifelse(y_bin < -model$tau, 0L, NA_integer_))
  y_lbl_con <- ifelse(y_bin < -model$tau, 1L, ifelse(y_bin >  model$tau, 0L, NA_integer_))
  sc_exp <- pred_bins
  sc_con <- -pred_bins

  # ---- accumulate for overall ROC ----
  keep_exp <- !is.na(y_lbl_exp)
  all_lbl_exp <- c(all_lbl_exp, y_lbl_exp[keep_exp])
  all_scr_exp <- c(all_scr_exp, sc_exp[keep_exp])

  keep_con <- !is.na(y_lbl_con)
  all_lbl_con <- c(all_lbl_con, y_lbl_con[keep_con])
  all_scr_con <- c(all_scr_con, sc_con[keep_con])

  # ROC plots (per-sample)
  roc_exp <- roc_plot_from(y_lbl_exp, sc_exp, "Expansion ROC (bin-level)")
  roc_con <- roc_plot_from(y_lbl_con, sc_con, "Contraction ROC (bin-level)")

  # AUCs for CSV
  get_auc <- function(lbl, sc) {
    ok <- !is.na(lbl) & !is.na(sc)
    if (length(unique(lbl[ok])) < 2) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(lbl[ok], sc[ok], quiet=TRUE, levels=c(0,1))))
  }
  auc_exp <- get_auc(y_lbl_exp, sc_exp)
  auc_con <- get_auc(y_lbl_con, sc_con)

  # PANEL: smoothed Δ, bin means, model predicted bin Δ (step), ROC
  dat <- dat %>% mutate(delta_s = roll_smooth(delta, k = smooth_k))
  bin_df <- tibble(bin = seq_along(mids), mid = mids, y_bin = y_bin, y_pred = pred_bins)

  plt_signal <- ggplot(dat, aes(x, delta_s)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = c(-tau, tau), linetype = "dotted", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(data = bin_df, aes(x = mid, y = y_bin), size = 1.2, alpha = 0.8, inherit.aes = FALSE) +
    geom_step(data = bin_df, aes(x = mid, y = y_pred), direction = "mid",
              linewidth = 0.8, inherit.aes = FALSE) +
    labs(title = bn,
         subtitle = paste0("Cluster ", cls, " (", cls_type, "); ",
                           "posteriors: ",
                           paste0("C", seq_along(probs), "=", scales::percent(probs, accuracy=0.1),
                                  collapse=", ")),
         x = "Repeat length", y = "Δ (Day42 − Day0)") +
    theme_minimal()

  panel <- plt_signal / (roc_exp | roc_con)

  print(panel)
  ggsave(file.path(out_dir, paste0(tools::file_path_sans_ext(bn), "_panel.png")),
         panel, width = 11, height = 8.5, units = "in", dpi = 150)

  metrics <- bind_rows(metrics, tibble(
    file = bn, method = method, K = model$K, tau = model$tau,
    assigned_cluster = cls, cluster_type = cls_type,
    auc_exp = auc_exp, auc_con = auc_con
  ))

  pred_bins_all <- bind_rows(pred_bins_all,
                             bin_df %>% mutate(file = bn, cluster = cls, cluster_type = cls_type))
}

# ===== Overall ROC page (minimalistic) =====
overall_exp_plot <- overall_roc_min_plot(all_lbl_exp, all_scr_exp, "Overall Expansion ROC (all test bins)")
overall_con_plot <- overall_roc_min_plot(all_lbl_con, all_scr_con, "Overall Contraction ROC (all test bins)")
overall_panel <- (overall_exp_plot | overall_con_plot)
print(overall_panel)

# Save metrics & per-bin predictions and close PDF containing ROC curve analysis
write_csv(metrics,       file.path(out_dir, "test_metrics.csv"))
write_csv(pred_bins_all, file.path(out_dir, "test_bin_predictions.csv"))
grDevices::dev.off()
message("Saved: ", pdf_fn) 
message("Saved: ", file.path(out_dir, "test_metrics.csv"))
message("Saved: ", file.path(out_dir, "test_bin_predictions.csv")) 
message("Saved per-file PNG panels to: ", out_dir)
