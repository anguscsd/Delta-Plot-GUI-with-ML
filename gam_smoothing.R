# ========================= USER INPUT =========================
folder      <- "/FOLDER/CONTAINING/DELTA/DATA"   # <- set this to your folder of CSVs
tau         <- 0.0075   # practical significance band around 0 for Δ (Day42-Day0)
min_len     <- 3        # minimum consecutive lengths to call a segment
k_spline    <- 20       # upper bound for GAM spline basis size
smooth_k    <- 5        # small pre-smoothing (rolling mean) window
cp_penalty  <- "MBIC"   # PELT penalty ("MBIC" is conservative)
iou_thr     <- 0.3      # IoU threshold for segment-level eval
max_length  <- 250      # only analyze repeat lengths <= this value
# ===============================================================

# ========================= SETUP ===============================
pkgs <- c("tidyverse","mgcv","gratia","changepoint","pROC",
          "zoo","patchwork","scales","readr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))

stopifnot(dir.exists(folder))
files <- list.files(folder, pattern="\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)

out_dir <- file.path(folder, "out")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ========================= HELPERS =============================
collapse_runs <- function(labels) {
  if (is.data.frame(labels)) labels <- dplyr::pull(labels, 1)
  if (is.list(labels))       labels <- unlist(labels, use.names = FALSE)
  if (is.factor(labels))     labels <- as.character(labels)
  labels <- as.vector(labels)
  r <- rle(labels)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  tibble(label = r$values, i1 = starts, i2 = ends, len = r$lengths)
}

label_from_ci <- function(lower, upper, tau) {
  ifelse(lower >  tau, "expansion",
  ifelse(upper < -tau, "contraction", "neutral"))
}

# ---- robust PELT labeling ----
cp_labels_vector <- function(y, tau, penalty="MBIC") {
  y2 <- zoo::rollmean(y, k = 3, align = "center", fill = "extend")
  cp <- changepoint::cpt.mean(y2, method = "PELT", penalty = penalty, class = TRUE)
  ends   <- c(changepoint::cpts(cp), length(y2))   # cpts are ends
  starts <- c(1, head(ends, -1) + 1)
  segs <- tibble(i1 = starts, i2 = ends)
  segs$seg_mean <- mapply(function(a,b) mean(y2[a:b], na.rm=TRUE), segs$i1, segs$i2)
  segs$label <- dplyr::case_when(
    segs$seg_mean >  tau ~ "expansion",
    segs$seg_mean < -tau ~ "contraction",
    TRUE ~ "neutral"
  )
  lab <- rep("neutral", length(y2))
  for (k in seq_len(nrow(segs))) lab[segs$i1[k]:segs$i2[k]] <- segs$label[k]
  lab
}

vec_to_segments <- function(x, v, labels, min_len = 3) {
  runs <- collapse_runs(labels) %>% dplyr::filter(label != "neutral", len >= min_len)
  if (nrow(runs) == 0) {
    return(tibble(label=character(), start_rl=numeric(), end_rl=numeric(),
                  len=integer(), mean_eff=numeric(), max_abs=numeric()))
  }
  tibble(
    label    = runs$label,
    start_rl = x[runs$i1],
    end_rl   = x[runs$i2],
    len      = runs$len,
    mean_eff = mapply(function(a,b) mean(v[a:b], na.rm=TRUE), runs$i1, runs$i2),
    max_abs  = mapply(function(a,b) max(abs(v[a:b]), na.rm=TRUE), runs$i1, runs$i2)
  )
}

truth_segments_from_y <- function(x, yf, tau, min_len) {
  lab <- ifelse(yf >  tau, "expansion",
         ifelse(yf < -tau, "contraction", "neutral"))
  vec_to_segments(x, v = yf, labels = lab, min_len = min_len)
}

interval_iou <- function(a1, a2, b1, b2) {
  inter <- pmax(0, pmin(a2, b2) - pmax(a1, b1))
  union <- pmax(a2, b2) - pmin(a1, b1)
  ifelse(union > 0, inter/union, 0)
}

segment_level_curves <- function(pred_segs, truth_segs, dir = c("expansion","contraction"),
                                 score_col = "mean_eff", iou_thr = 0.3) {
  dir <- match.arg(dir)
  P <- pred_segs %>% dplyr::filter(label == dir) %>%
       dplyr::mutate(score = abs(.data[[score_col]]))
  T <- truth_segs %>% dplyr::filter(label == dir)

  if (nrow(P) == 0 || nrow(T) == 0) {
    return(list(
      roc = tibble(fpr = c(0,1), tpr = c(0,1)),
      pr  = tibble(recall = c(0,1), precision = c(1, nrow(T)/(nrow(T)+1e-9))),
      auc_roc = NA_real_, auc_pr = NA_real_,
      n_pred = nrow(P), n_truth = nrow(T)
    ))
  }

  best_iou <- sapply(seq_len(nrow(P)), function(i) {
    ious <- interval_iou(P$start_rl[i], P$end_rl[i], T$start_rl, T$end_rl)
    max(ious, na.rm = TRUE)
  })
  P$best_iou <- best_iou

  thrs <- sort(unique(P$score), decreasing = TRUE)
  curve_df <- tibble(threshold = thrs, TP = NA_integer_, FP = NA_integer_, FN = NA_integer_)

  for (k in seq_along(thrs)) {
    thr <- thrs[k]
    S <- P %>% dplyr::filter(score >= thr)

    is_TP <- if (nrow(S) > 0) S$best_iou >= iou_thr else logical(0)
    TP <- sum(is_TP)
    FP <- nrow(S) - TP

    matched_T <- rep(FALSE, nrow(T))
    if (nrow(S) > 0) {
      for (i in which(is_TP)) {
        ious <- interval_iou(S$start_rl[i], S$end_rl[i], T$start_rl, T$end_rl)
        matched_T <- matched_T | (ious >= iou_thr)
      }
    }
    FN <- sum(!matched_T)

    curve_df$TP[k] <- TP; curve_df$FP[k] <- FP; curve_df$FN[k] <- FN
  }

  curve_df <- curve_df %>%
    dplyr::mutate(
      precision = TP / pmax(TP + FP, 1e-9),
      recall    = TP / pmax(TP + FN, 1e-9),
      tpr       = recall,
      fpr       = FP / pmax(FP + (nrow(T) - TP), 1e-9)
    ) %>%
    dplyr::arrange(dplyr::desc(threshold))

  auc_trap <- function(x, y) {
    o <- order(x, decreasing = FALSE)
    if (length(o) < 2) return(NA_real_)
    sum(diff(x[o]) * zoo::rollmean(y[o], 2))
  }
  auc_roc <- auc_trap(curve_df$fpr, curve_df$tpr)
  auc_pr  <- auc_trap(curve_df$recall, curve_df$precision)

  list(
    roc = curve_df %>% dplyr::select(fpr, tpr),
    pr  = curve_df %>% dplyr::select(recall, precision),
    auc_roc = auc_roc,
    auc_pr  = auc_pr,
    n_pred = nrow(P),
    n_truth = nrow(T)
  )
}

plot_pr <- function(df, auc, title_txt, note = NULL) {
  if (any(is.na(auc))) {
    return(
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "PR unavailable", size = 4) +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = title_txt, x = "Recall", y = "Precision") + theme_minimal()
    )
  }
  ann <- sprintf("AUPR = %.3f%s", auc, if (!is.null(note)) paste0("\n", note) else "")
  ggplot(df, aes(recall, precision)) +
    geom_path() +
    annotate("text", x = 0.65, y = 0.1, hjust = 0, label = ann) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    labs(title = title_txt, x = "Recall", y = "Precision") +
    theme_minimal()
}

plot_roc <- function(df, auc, title_txt, note = NULL) {
  if (any(is.na(auc))) {
    return(
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "ROC unavailable", size = 4) +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = title_txt, x = "False Positive Rate", y = "True Positive Rate") +
        theme_minimal()
    )
  }
  ann <- sprintf("AUC = %.3f%s", auc, if (!is.null(note)) paste0("\n", note) else "")
  ggplot(df, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_path() +
    annotate("text", x = 0.65, y = 0.1, hjust = 0, label = ann) +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(title = title_txt, x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal()
}

# ========================= CORE ANALYSIS (SINGLE FILE) =====================
analyze_one_file <- function(csv_path,
                             tau = 0.003,
                             min_len = 5,
                             k_spline = 60,
                             smooth_k = 7,
                             cp_penalty = "MBIC",
                             iou_thr = 0.3,
                             max_length = 200) {

  df <- read.csv(csv_path, check.names = FALSE)
  stopifnot("delta" %in% names(df))
  length_col <- names(df)[1]

  # ---- NEW: sort then hard-trim to max_length BEFORE any processing ----
  df <- df %>% dplyr::arrange(.data[[length_col]]) %>%
       dplyr::filter(.data[[length_col]] <= max_length)

  # If everything got trimmed, bail gracefully
  if (nrow(df) < 3) {
    empty_panel <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(basename(csv_path), "\n(no data ≤ ", max_length, ")"),
               size = 4) +
      xlim(0,1) + ylim(0,1) + theme_void()
    return(list(
      plot = empty_panel,
      seg_gam = tibble(), seg_cp = tibble(), seg_consensus = tibble(),
      seg_truth = tibble(),
      metrics = tibble(
        file = basename(csv_path), n_points = 0,
        tau = tau, min_len = min_len, k_spline = NA_integer_, iou_thr = iou_thr,
        n_truth_exp = 0, n_pred_exp = 0, n_truth_con = 0, n_pred_con = 0,
        aupr_exp = NA_real_, auroc_exp = NA_real_, aupr_con = NA_real_, auroc_con = NA_real_,
        summary = paste0("no data ≤ ", max_length)
      )
    ))
  }

  x  <- as.numeric(df[[length_col]])
  y  <- as.numeric(df$delta)

  # pre-smooth (gentle)
  ys <- zoo::rollmean(y, k = smooth_k, align = "center", fill = NA)
  yf <- ifelse(is.na(ys), y, ys)

  # GAM fit
  n <- length(x)
  k_use <- max(10, min(k_spline, floor(n/4)))
  gam_fit <- mgcv::gam(yf ~ s(x, k = k_use), method = "REML")

  # predictions & simultaneous CIs at observed x
  pred <- predict(gam_fit, newdata = tibble(x=x), se.fit = TRUE, type = "response")
  crit <- tryCatch({
    cs <- gratia::crit_smooth(gam_fit, coverage = 0.95)
    if ("smooth" %in% names(cs)) cs$crit[match("s(x)", cs$smooth)] else cs$crit[1]
  }, error = function(e) NA_real_)
  if (is.na(crit) || !is.finite(crit)) crit <- qnorm(0.975)

  lower <- pred$fit - crit*pred$se.fit
  upper <- pred$fit + crit*pred$se.fit
  lab_gam <- label_from_ci(lower, upper, tau)

  # Change-point labels
  lab_cp <- cp_labels_vector(yf, tau = tau, penalty = cp_penalty)

  # Consensus per index
  lab_consensus <- ifelse(lab_gam == "expansion" & lab_cp == "expansion", "expansion",
                    ifelse(lab_gam == "contraction" & lab_cp == "contraction", "contraction",
                           "neutral"))

  # Segment tables
  seg_gam       <- vec_to_segments(x, pred$fit, lab_gam,       min_len = min_len)
  seg_cp        <- vec_to_segments(x, yf,       lab_cp,        min_len = min_len)
  seg_consensus <- vec_to_segments(x, pred$fit, lab_consensus, min_len = min_len)

  # Truth segments from smoothed raw Δ
  seg_truth <- truth_segments_from_y(x, yf, tau = tau, min_len = min_len)

  # One-line summary
  summary_text <- if (nrow(seg_consensus) == 0) {
    "no robust change segments"
  } else {
    paste(
      seg_consensus %>%
        dplyr::group_by(label) %>%
        dplyr::summarise(ranges = paste0(min(start_rl), "–", max(end_rl), collapse="; "), .groups="drop") %>%
        dplyr::mutate(txt = paste0(label, " between lengths ", ranges)) %>%
        dplyr::pull(txt),
      collapse = ", "
    )
  }

  # -------------- SEGMENT-LEVEL PR & ROC (for metrics only) ---------------
  curves_exp <- segment_level_curves(seg_consensus, seg_truth, dir = "expansion",
                                     score_col = "mean_eff", iou_thr = iou_thr)
  curves_con <- segment_level_curves(seg_consensus, seg_truth, dir = "contraction",
                                     score_col = "mean_eff", iou_thr = iou_thr)

  # -------------- Plot: signal + consensus band (NO PR/ROC in output) -----
  rng_y <- range(c(y, lower, upper), na.rm = TRUE)
  y_span <- diff(rng_y); if (!is.finite(y_span) || y_span == 0) y_span <- 1
  y_bottom <- rng_y[1] - 0.12*y_span
  band_h   <- 0.06*y_span

  cons_band <- if (nrow(seg_consensus) > 0) {
    seg_consensus %>% dplyr::transmute(
      xmin = start_rl, xmax = end_rl,
      ymin = y_bottom, ymax = y_bottom + band_h,
      label = label
    )
  } else tibble(xmin=numeric(), xmax=numeric(), ymin=numeric(), ymax=numeric(), label=character())

  plt_signal <- ggplot(tibble(x, y, yf, fit = pred$fit, lower, upper, label = lab_gam), aes(x, fit)) +
    { if (nrow(cons_band) > 0)
      geom_rect(data = cons_band,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
                inherit.aes = FALSE, alpha = 0.35) } +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = c(-tau, tau), linetype = "dotted", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(x = x, y = y), size = 0.6, alpha = 0.7, inherit.aes = FALSE) +
    scale_fill_manual(values = c(expansion = "#1f77b4", contraction = "#d62728")) +
    labs(title = basename(csv_path),
         subtitle = summary_text,
         x = "Repeat length",
         y = "Δ (Day42 − Day0)",
         fill = "Consensus") +
    coord_cartesian(ylim = c(y_bottom, rng_y[2])) +
    theme_classic()

  # Only the signal panel goes to outputs
  panel <- plt_signal

  list(
    plot          = panel,
    seg_gam       = seg_gam %>% dplyr::mutate(file = basename(csv_path), .before=1),
    seg_cp        = seg_cp  %>% dplyr::mutate(file = basename(csv_path), .before=1),
    seg_consensus = seg_consensus %>% dplyr::mutate(file = basename(csv_path), .before=1),
    seg_truth     = seg_truth %>% dplyr::mutate(file = basename(csv_path), .before=1),
    metrics = tibble(
      file = basename(csv_path),
      n_points = length(x),
      tau = tau, min_len = min_len, k_spline = k_use,
      iou_thr = iou_thr, max_length = max_length,
      n_truth_exp = curves_exp$n_truth, n_pred_exp = curves_exp$n_pred,
      n_truth_con = curves_con$n_truth, n_pred_con = curves_con$n_pred,
      aupr_exp = curves_exp$auc_pr, auroc_exp = curves_exp$auc_roc,
      aupr_con = curves_con$auc_pr, auroc_con = curves_con$auc_roc,
      summary = summary_text
    )
  )
}

# ========================= BATCH OVER FOLDER ==============================
all_metrics <- tibble()
all_segs    <- tibble()
plot_list   <- list()  #
pdf_fn <- file.path(out_dir, "GAM_report.pdf") # change name or location of output report if required
grDevices::pdf(pdf_fn, width = 11, height = 8.5)
on.exit(grDevices::dev.off(), add = TRUE)

message("Processing ", length(files), " file(s)...")
for (f in files) {
  res <- analyze_one_file(f, tau = tau, min_len = min_len,
                          k_spline = k_spline, smooth_k = smooth_k,
                          cp_penalty = cp_penalty, iou_thr = iou_thr,
                          max_length = max_length)

  # Save reference to the plot variable
  plot_var_name <- paste0("plot_", tools::file_path_sans_ext(basename(f)))
  plot_list[[plot_var_name]] <- res$plot

  # Write to PDF
  print(res$plot)

  # Save per-file PNG
  png_fn <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), "_panel.png"))
  ggsave(png_fn, res$plot, width = 11, height = 8.5, units = "in", dpi = 150)

  # Collect metrics and segments
  all_metrics <- dplyr::bind_rows(all_metrics, res$metrics)
  if (nrow(res$seg_consensus) > 0) {
    all_segs <- dplyr::bind_rows(all_segs, res$seg_consensus %>% dplyr::mutate(source="consensus"))
  } else {
    all_segs <- dplyr::bind_rows(all_segs, tibble(
      file = basename(f), label="none", start_rl=NA_real_, end_rl=NA_real_,
      len=NA_integer_, mean_eff=NA_real_, max_abs=NA_real_, source="consensus"
    ))
  }
  all_segs <- dplyr::bind_rows(
    all_segs,
    res$seg_gam %>% dplyr::mutate(source = "GAM"),
    res$seg_cp  %>% dplyr::mutate(source = "PELT")
  )
}

# ========================= SAVE CSV SUMMARIES =============================
metrics_csv <- file.path(out_dir, "metrics.csv")
segs_csv    <- file.path(out_dir, "segments.csv")
readr::write_csv(all_metrics, metrics_csv)
readr::write_csv(all_segs,    segs_csv)

# ========================= PRINT PLOT LIST =============================
# Useful for downstream analysis or figure generation of specific datasets/ visualisations
message("\nPlots generated:")
for (i in seq_along(plot_list)) {
  cat(sprintf("%-40s -> %s\n",
              names(plot_list)[i],
              names(plot_list)[i]))
}

message("\nSaved: ", pdf_fn)
message("Saved: ", metrics_csv)
message("Saved: ", segs_csv)
