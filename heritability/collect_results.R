#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x) || identical(x, "")) y else x
parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Invalid argument near: ", x[[i]])
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
sources <- strsplit(args[["sources"]] %||% "pangenie", ",", fixed = TRUE)[[1]]
out_base <- args[["out-base"]] %||% stop("Missing --out-base")
summary_dir <- args[["summary-dir"]] %||% file.path(out_base, "summary")
plot_dir <- args[["plot-dir"]] %||% file.path(out_base, "plots")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

write_both <- function(x, stem) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, paste0(stem, ".tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(x, paste0(stem, ".csv"), quote = TRUE, na = "NA")
}

task_list <- list(); status_list <- list()
for (source in sources) {
  task_file <- file.path(out_base, source, "inputs", paste0(source, ".reml_tasks.tsv"))
  if (!file.exists(task_file)) stop("Missing task manifest: ", task_file)
  t <- fread(task_file)
  if (nrow(t) != 34L) warning(source, " task manifest contains ", nrow(t), " rows; expected 34")
  task_list[[source]] <- t
  status_file <- file.path(out_base, source, "status", paste0(source, ".reml_status.tsv"))
  if (file.exists(status_file)) status_list[[source]] <- fread(status_file, fill = TRUE)
}
tasks <- rbindlist(task_list, fill = TRUE)
if (anyDuplicated(tasks[, .(source, task_id)])) stop("Duplicated source/task_id in task manifests")
statuses <- if (length(status_list)) rbindlist(status_list, fill = TRUE) else data.table()
if (nrow(statuses) && anyDuplicated(statuses[, .(source, task_id)])) {
  setorder(statuses, source, task_id, finished_at)
  statuses <- statuses[, .SD[.N], by = .(source, task_id)]
}
status_cols <- c("source", "task_id", "status", "exit_code", "error_reason", "finished_at")
if (!nrow(statuses)) statuses <- data.table(source = character(), task_id = character(), status = character(),
  exit_code = integer(), error_reason = character(), finished_at = character())
statuses <- statuses[, intersect(status_cols, names(statuses)), with = FALSE]
task_status <- merge(tasks, statuses, by = c("source", "task_id"), all.x = TRUE)
task_status[is.na(status), `:=`(status = "not_run", error_reason = "No status record", exit_code = NA_integer_)]
task_status[, hsq_file := paste0(out_prefix, ".hsq")]

map_component <- function(source_row, components) {
  comps <- strsplit(components, ",", fixed = TRUE)[[1]]
  if (grepl("^V\\(G\\)(/Vp|/Vp_L)$", source_row)) return(comps[[1]])
  m <- regexec("^V\\(G([0-9]+)\\)", source_row)
  hit <- regmatches(source_row, m)[[1]]
  if (length(hit) >= 2L) {
    idx <- as.integer(hit[[2]])
    if (!is.na(idx) && idx <= length(comps)) return(comps[[idx]])
  }
  NA_character_
}

read_one_hsq <- function(task) {
  file <- task$hsq_file
  if (!file.exists(file) || file.info(file)$size == 0) {
    return(data.table(Source = NA_character_, Variance = NA_real_, SE = NA_real_, component = NA_character_))
  }
  x <- tryCatch(fread(file, fill = TRUE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(data.table(Source = NA_character_, Variance = NA_real_, SE = NA_real_, component = NA_character_))
  if (!"Source" %in% names(x)) setnames(x, 1, "Source")
  if (!"Variance" %in% names(x) && ncol(x) >= 2L) setnames(x, 2, "Variance")
  if (!"SE" %in% names(x) && ncol(x) >= 3L) setnames(x, 3, "SE")
  if (!"Variance" %in% names(x)) x[, Variance := NA_real_]
  if (!"SE" %in% names(x)) x[, SE := NA_real_]
  x[, `:=`(Source = as.character(Source), Variance = suppressWarnings(as.numeric(Variance)), SE = suppressWarnings(as.numeric(SE)))]
  x[, component := vapply(Source, map_component, character(1), task$components)]
  x
}

long_list <- vector("list", nrow(task_status))
for (i in seq_len(nrow(task_status))) {
  task <- task_status[i]
  h <- read_one_hsq(task)
  meta <- task[, .(source, task_id, model, pc_n, adjustment, components, status, exit_code,
    error_reason, input_type, input_path, qcovar, out_prefix, log_file, hsq_file)]
  long_list[[i]] <- cbind(meta[rep(1L, nrow(h))], h)
}
long <- rbindlist(long_list, fill = TRUE)
write_both(long, file.path(summary_dir, "heritability_all_models.long"))

first_value <- function(h, pattern, field = "Variance", fixed = FALSE) {
  idx <- if (fixed) which(h$Source %chin% pattern) else grep(pattern, h$Source)
  if (!length(idx)) return(NA_real_)
  as.numeric(h[[field]][idx[[1]]])
}

make_total <- function(task, h) {
  comps <- strsplit(task$components, ",", fixed = TRUE)[[1]]
  is_joint <- length(comps) > 1L
  obs_total <- if (is_joint) first_value(h, "Sum of V(G)/Vp", fixed = TRUE) else first_value(h, "V(G)/Vp", fixed = TRUE)
  obs_se <- if (is_joint) first_value(h, "Sum of V(G)/Vp", "SE", TRUE) else first_value(h, "V(G)/Vp", "SE", TRUE)
  liab_total <- if (is_joint) first_value(h, c("Sum of V(G)_L/Vp", "Sum of V(G)/Vp_L"), fixed = TRUE) else first_value(h, "V(G)/Vp_L", fixed = TRUE)
  liab_se <- if (is_joint) first_value(h, c("Sum of V(G)_L/Vp", "Sum of V(G)/Vp_L"), "SE", TRUE) else first_value(h, "V(G)/Vp_L", "SE", TRUE)
  obs_rows <- h[grepl("^V\\(G[0-9]*\\)/Vp$", Source)]
  liab_rows <- h[grepl("^V\\(G[0-9]*\\)/Vp_L$", Source)]
  if (!is.finite(obs_total) && nrow(obs_rows)) obs_total <- sum(obs_rows$Variance, na.rm = TRUE)
  if (!is.finite(liab_total) && nrow(liab_rows)) liab_total <- sum(liab_rows$Variance, na.rm = TRUE)
  out <- data.table(
    source = task$source, task_id = task$task_id, model = task$model,
    pc_n = task$pc_n, adjustment = task$adjustment, components = task$components,
    status = task$status, exit_code = task$exit_code, error_reason = task$error_reason,
    n = first_value(h, "n", fixed = TRUE), p_value = first_value(h, "Pval", fixed = TRUE),
    h2_observed_total = obs_total, se_observed_total = obs_se,
    h2_liability_total = liab_total, se_liability_total = liab_se,
    h2_snv_observed = NA_real_, h2_sv_observed = NA_real_, h2_tr_observed = NA_real_,
    h2_snv_liability = NA_real_, h2_sv_liability = NA_real_, h2_tr_liability = NA_real_,
    log_file = task$log_file, hsq_file = task$hsq_file
  )
  for (comp in c("SNV_INDEL", "SV", "TR")) {
    ov <- obs_rows[component == comp, Variance][1]
    lv <- liab_rows[component == comp, Variance][1]
    short <- c(SNV_INDEL = "snv", SV = "sv", TR = "tr")[[comp]]
    if (length(ov)) set(out, j = paste0("h2_", short, "_observed"), value = as.numeric(ov))
    if (length(lv)) set(out, j = paste0("h2_", short, "_liability"), value = as.numeric(lv))
  }
  out
}

totals_list <- vector("list", nrow(task_status))
for (i in seq_len(nrow(task_status))) totals_list[[i]] <- make_total(task_status[i], read_one_hsq(task_status[i]))
totals <- rbindlist(totals_list, fill = TRUE)
write_both(totals, file.path(summary_dir, "heritability_model_totals"))

focus_map <- c(
  "SNV_INDEL.marginal" = "SNV",
  "SNV_INDEL_SV.joint" = "SNV+SV",
  "SNV_INDEL_TR.joint" = "SNV+TR",
  "SNV_INDEL_SV_TR.joint" = "SNV+SV+TR"
)
focus <- totals[model %chin% names(focus_map)]
focus[, focus_model := unname(focus_map[model])]
focus[, focus_model := factor(focus_model, levels = c("SNV", "SNV+SV", "SNV+TR", "SNV+SV+TR"))]
setcolorder(focus, c("source", "focus_model", setdiff(names(focus), c("source", "focus_model"))))
write_both(focus, file.path(summary_dir, "heritability_focus_snv_models.long"))

make_focus_wide <- function(x, selection_name) {
  if (!nrow(x)) return(data.table())
  value_cols <- c("h2_observed_total", "se_observed_total", "h2_liability_total", "se_liability_total", "status", "error_reason", "n")
  w <- dcast(x, source + pc_n ~ focus_model, value.var = value_cols)
  names(w) <- gsub("\\+", "_plus_", tolower(names(w)))
  names(w) <- gsub("[^a-z0-9_]+", "_", names(w))
  obs <- function(model) paste0("h2_observed_total_", model)
  liab <- function(model) paste0("h2_liability_total_", model)
  if (all(c(obs("snv"), obs("snv_plus_sv")) %in% names(w))) w[, delta_observed_add_sv := get(obs("snv_plus_sv")) - get(obs("snv"))]
  if (all(c(obs("snv"), obs("snv_plus_tr")) %in% names(w))) w[, delta_observed_add_tr := get(obs("snv_plus_tr")) - get(obs("snv"))]
  if (all(c(obs("snv"), obs("snv_plus_sv_plus_tr")) %in% names(w))) w[, delta_observed_add_sv_tr := get(obs("snv_plus_sv_plus_tr")) - get(obs("snv"))]
  if (all(c(liab("snv"), liab("snv_plus_sv")) %in% names(w))) w[, delta_liability_add_sv := get(liab("snv_plus_sv")) - get(liab("snv"))]
  if (all(c(liab("snv"), liab("snv_plus_tr")) %in% names(w))) w[, delta_liability_add_tr := get(liab("snv_plus_tr")) - get(liab("snv"))]
  if (all(c(liab("snv"), liab("snv_plus_sv_plus_tr")) %in% names(w))) w[, delta_liability_add_sv_tr := get(liab("snv_plus_sv_plus_tr")) - get(liab("snv"))]
  w[, comparison := selection_name]
  setcolorder(w, c("source", "pc_n", "comparison", setdiff(names(w), c("source", "pc_n", "comparison"))))
  w
}

baseline <- focus[adjustment == "snv_pc"]
baseline_wide <- make_focus_wide(baseline, "same SNV PCs in all four models")
write_both(baseline_wide, file.path(summary_dir, "heritability_focus_snv_pc_only.wide"))

fully <- focus[
  (focus_model == "SNV" & adjustment == "snv_pc") |
  (focus_model == "SNV+SV" & adjustment == "snv_sv_pc") |
  (focus_model == "SNV+TR" & adjustment == "snv_tr_pc") |
  (focus_model == "SNV+SV+TR" & adjustment == "snv_sv_tr_pc")
]
fully_wide <- make_focus_wide(fully, "component-matched PCs; covariates differ between models")
write_both(fully_wide, file.path(summary_dir, "heritability_focus_fully_adjusted.wide"))

for (src in sources) {
  source_dir <- file.path(summary_dir, "by_source", src)
  write_both(totals[source == src], file.path(source_dir, paste0(src, ".all_model_totals")))
  write_both(focus[source == src], file.path(source_dir, paste0(src, ".focus_snv_models.long")))
  write_both(baseline_wide[source == src], file.path(source_dir, paste0(src, ".focus_snv_pc_only.wide")))
  write_both(fully_wide[source == src], file.path(source_dir, paste0(src, ".focus_fully_adjusted.wide")))
}

errors <- totals[!status %chin% c("success", "reused")]
write_both(errors, file.path(summary_dir, "heritability_errors"))

plot_focus <- focus[status %chin% c("success", "reused", "warning") & is.finite(h2_liability_total)]
if (nrow(plot_focus)) {
  p <- ggplot(plot_focus, aes(x = focus_model, y = h2_liability_total, colour = adjustment, group = adjustment)) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_errorbar(aes(ymin = h2_liability_total - se_liability_total,
      ymax = h2_liability_total + se_liability_total), width = 0.16, na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.35), size = 2.1) +
    facet_grid(source ~ pc_n, scales = "free_y", labeller = label_both) +
    labs(title = "Figure7 heritability comparison", subtitle = "Liability-scale estimates; error bars are GCTA SE",
      x = NULL, y = "Liability-scale heritability", colour = "PC adjustment") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(plot_dir, "heritability_focus_comparison.png"), p, width = 12, height = 12, dpi = 180)
  ggsave(file.path(plot_dir, "heritability_focus_comparison.pdf"), p, width = 12, height = 12)

  pb <- ggplot(baseline[status %chin% c("success", "reused", "warning")],
      aes(x = focus_model, y = h2_liability_total, colour = source, group = source)) +
    geom_errorbar(aes(ymin = h2_liability_total - se_liability_total,
      ymax = h2_liability_total + se_liability_total), width = 0.14, na.rm = TRUE) +
    geom_point(size = 2.2) + facet_wrap(~pc_n, labeller = label_both) +
    labs(title = "Heritability comparison with the same SNV PC adjustment",
      x = NULL, y = "Liability-scale heritability", colour = "Source") +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(plot_dir, "heritability_focus_snv_pc_only.png"), pb, width = 10, height = 6.5, dpi = 180)
  ggsave(file.path(plot_dir, "heritability_focus_snv_pc_only.pdf"), pb, width = 10, height = 6.5)
}

liab_component_rows <- long[
  status %chin% c("success", "reused", "warning") &
    is.finite(Variance) &
    (
      (grepl("\\.marginal$", model) & Source == "V(G)/Vp_L") |
        (grepl("\\.joint$", model) & grepl("^V\\(G[0-9]+\\)/Vp_L$", Source))
    )
]
if (nrow(liab_component_rows)) {
  model_labels <- c(
    SNV_INDEL.marginal = "SNV",
    SV.marginal = "SV",
    TR.marginal = "TR",
    SNV_INDEL_SV.joint = "SNV+SV",
    SNV_INDEL_TR.joint = "SNV+TR",
    SV_TR.joint = "SV+TR",
    SNV_INDEL_SV_TR.joint = "SNV+SV+TR"
  )
  adjustment_labels <- c(
    snv_pc = "SNV PC",
    sv_pc = "SV PC",
    snv_sv_pc = "SNV+SV PC",
    tr_pc = "TR PC",
    snv_tr_pc = "SNV+TR PC",
    sv_tr_pc = "SV+TR PC",
    snv_sv_tr_pc = "SNV+SV+TR PC"
  )
  liab_component_rows[, `:=`(
    model_label = unname(model_labels[model]),
    adjustment_label = unname(adjustment_labels[adjustment]),
    component_label = fifelse(component == "SNV_INDEL", "SNV_INDEL", as.character(component)),
    h2_liability_component = Variance,
    se_liability_component = SE
  )]
  liab_component_rows[is.na(model_label), model_label := model]
  liab_component_rows[is.na(adjustment_label), adjustment_label := adjustment]
  liab_component_rows[is.na(component_label), component_label := "Unknown"]
  liab_component_rows[, model_label := factor(model_label,
    levels = c("SNV", "SV", "TR", "SNV+SV", "SNV+TR", "SV+TR", "SNV+SV+TR"))]
  liab_component_rows[, component_label := factor(component_label, levels = c("SNV_INDEL", "SV", "TR", "Unknown"))]
  liab_component_rows[, adjustment_label := factor(adjustment_label,
    levels = c("SNV PC", "SV PC", "SNV+SV PC", "TR PC", "SNV+TR PC", "SV+TR PC", "SNV+SV+TR PC"))]
  setorder(liab_component_rows, source, pc_n, adjustment_label, model_label, component_label)
  write_both(liab_component_rows, file.path(summary_dir, "heritability_liability_components_for_plot"))

  component_colours <- c(SNV_INDEL = "#4C78A8", SV = "#F58518", TR = "#54A24B", Unknown = "grey65")
  make_liab_component_plot <- function(pc_value) {
    x <- liab_component_rows[pc_n == pc_value]
    ggplot(x, aes(x = model_label, y = h2_liability_component, fill = component_label)) +
      geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
      geom_col(width = 0.72, colour = "grey30", linewidth = 0.15, na.rm = TRUE) +
      facet_grid(source ~ adjustment_label, scales = "free_x", space = "free_x", drop = TRUE) +
      scale_fill_manual(values = component_colours, drop = FALSE) +
      labs(
        title = sprintf("Liability-scale heritability components (PC%d)", pc_value),
        subtitle = "Marginal models use V(G)/Vp_L; joint models stack V(G1)/Vp_L, V(G2)/Vp_L, and V(G3)/Vp_L by component",
        x = NULL, y = "Liability-scale heritability", fill = "GRM component"
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        panel.spacing.x = unit(0.55, "lines"),
        panel.spacing.y = unit(0.7, "lines"),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "bottom"
      )
  }

  liab_component_plots <- list()
  for (pc_value in sort(unique(liab_component_rows$pc_n))) {
    pcomp <- make_liab_component_plot(pc_value)
    liab_component_plots[[paste0("pc", pc_value)]] <- pcomp
    ggsave(file.path(plot_dir, sprintf("heritability_liability_components_pc%d.png", pc_value)),
      pcomp, width = 22, height = 12, dpi = 180, limitsize = FALSE)
    ggsave(file.path(plot_dir, sprintf("heritability_liability_components_pc%d.pdf", pc_value)),
      pcomp, width = 22, height = 12, limitsize = FALSE)
  }
  pdf(file.path(plot_dir, "heritability_liability_components_all.pdf"), width = 22, height = 12, onefile = TRUE)
  for (pcomp in liab_component_plots) print(pcomp)
  dev.off()

  selected_specs <- data.table(
    model = c("SNV_INDEL.marginal", "SV.marginal", "TR.marginal",
      "SNV_INDEL_SV.joint", "SNV_INDEL_SV_TR.joint"),
    adjustment = c("snv_pc", "snv_sv_pc", "snv_tr_pc", "snv_sv_pc", "snv_sv_tr_pc"),
    selected_type = c("SNV", "SV", "TR", "SNV+SV", "SNV+SV+TR"),
    selected_pc = c("SNV PC", "SNV+SV PC", "SNV+TR PC", "SNV+SV PC", "SNV+SV+TR PC")
  )
  liab_selected <- merge(liab_component_rows, selected_specs,
    by = c("model", "adjustment"), all = FALSE, allow.cartesian = TRUE)
  if (nrow(liab_selected)) {
    liab_selected[, selected_type := factor(selected_type,
      levels = c("SNV", "SV", "TR", "SNV+SV", "SNV+SV+TR"))]
    liab_selected[, selected_pc := factor(selected_pc,
      levels = c("SNV PC", "SNV+SV PC", "SNV+TR PC", "SNV+SV+TR PC"))]
    liab_selected[, plot_label := factor(paste0(as.character(selected_type), "\n", as.character(selected_pc)),
      levels = paste0(c("SNV", "SV", "TR", "SNV+SV", "SNV+SV+TR"), "\n",
        c("SNV PC", "SNV+SV PC", "SNV+TR PC", "SNV+SV PC", "SNV+SV+TR PC")))]
    setorder(liab_selected, source, pc_n, selected_type, component_label)
    write_both(liab_selected, file.path(summary_dir, "heritability_liability_selected5_for_plot"))

    psel <- ggplot(liab_selected, aes(x = plot_label, y = h2_liability_component, fill = component_label)) +
      geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
      geom_col(width = 0.72, colour = "grey30", linewidth = 0.15, na.rm = TRUE) +
      facet_grid(source ~ pc_n, scales = "free_y", labeller = label_both) +
      scale_fill_manual(values = component_colours, drop = FALSE) +
      labs(
        title = "Selected liability-scale heritability components",
        subtitle = "Selected PC strategy: SNV=SNV PC; SV=SNV+SV PC; TR=SNV+TR PC; SNV+SV=SNV+SV PC; SNV+SV+TR=SNV+SV+TR PC",
        x = NULL, y = "Liability-scale heritability", fill = "GRM component"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        panel.spacing.x = unit(0.7, "lines"),
        panel.spacing.y = unit(0.75, "lines"),
        legend.position = "bottom"
      )
    ggsave(file.path(plot_dir, "heritability_liability_selected5_components.png"),
      psel, width = 13, height = 11, dpi = 180, limitsize = FALSE)
    ggsave(file.path(plot_dir, "heritability_liability_selected5_components.pdf"),
      psel, width = 13, height = 11, limitsize = FALSE)
  }

  overview_model_labels <- c(
    SNV_INDEL.marginal = "SNV",
    SV.marginal = "SV",
    TR.marginal = "TR",
    SNV_INDEL_SV.joint = "SNV+SV",
    SNV_INDEL_TR.joint = "SNV+TR",
    SV_TR.joint = "TR+SV",
    SNV_INDEL_SV_TR.joint = "SNV+SV+TR"
  )
  overview_adjustment_labels <- c(
    snv_pc = "SNV PC",
    sv_pc = "SV PC",
    tr_pc = "TR PC",
    snv_sv_pc = "SNV+SV PC",
    snv_tr_pc = "SNV+TR PC",
    sv_tr_pc = "TR+SV PC",
    snv_sv_tr_pc = "SNV+SV+TR PC"
  )
  overview_model_levels <- c("SNV", "SV", "TR", "SNV+SV", "SNV+TR", "TR+SV", "SNV+SV+TR")
  overview_adjustment_levels <- c("SNV PC", "SV PC", "TR PC", "SNV+SV PC",
    "SNV+TR PC", "TR+SV PC", "SNV+SV+TR PC")
  liab_overview <- copy(liab_component_rows)
  liab_overview[, `:=`(
    overview_model = unname(overview_model_labels[model]),
    overview_pc = unname(overview_adjustment_labels[adjustment])
  )]
  liab_overview <- liab_overview[!is.na(overview_model) & !is.na(overview_pc)]
  if (nrow(liab_overview)) {
    liab_overview[, `:=`(
      overview_model = factor(overview_model, levels = overview_model_levels),
      overview_pc = factor(overview_pc, levels = overview_adjustment_levels)
    )]
    overview_combos <- unique(liab_overview[, .(overview_model, overview_pc)])
    setorder(overview_combos, overview_model, overview_pc)
    overview_combos[, overview_x := factor(
      paste0(as.character(overview_model), "\n", as.character(overview_pc)),
      levels = paste0(as.character(overview_model), "\n", as.character(overview_pc))
    )]
    liab_overview <- merge(liab_overview, overview_combos,
      by = c("overview_model", "overview_pc"), all.x = TRUE, sort = FALSE)
    setorder(liab_overview, pc_n, source, overview_model, overview_pc, component_label)
    write_both(liab_overview, file.path(summary_dir, "heritability_liability_overview_unfaceted_for_plot"))

    make_overview_plot <- function(pc_value) {
      x <- liab_overview[pc_n == pc_value]
      p <- ggplot(x, aes(x = overview_x, y = h2_liability_component, fill = component_label)) +
        geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
        geom_col(width = 0.72, colour = "grey30", linewidth = 0.15, na.rm = TRUE) +
        scale_fill_manual(values = component_colours, drop = FALSE) +
        labs(
          title = sprintf("Liability-scale heritability overview (PC%d)", pc_value),
          subtitle = if (uniqueN(x$source) > 1L)
            "X-axis ordered by model then PC adjustment; panels are sources; bars are stacked by GRM component"
          else
            "X-axis ordered by model then PC adjustment; bars are stacked by GRM component",
          x = NULL, y = "Liability-scale heritability", fill = "GRM component"
        ) +
        theme_bw(base_size = 10.5) +
        theme(
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom"
        )
      if (uniqueN(x$source) > 1L) p <- p + facet_grid(source ~ ., scales = "free_y")
      p
    }
    overview_plots <- list()
    overview_width <- 15
    overview_height <- if (uniqueN(liab_overview$source) > 1L) 12 else 8.5
    for (pc_value in sort(unique(liab_overview$pc_n))) {
      poverview <- make_overview_plot(pc_value)
      overview_plots[[paste0("pc", pc_value)]] <- poverview
      ggsave(file.path(plot_dir, sprintf("heritability_liability_overview_unfaceted_pc%d.png", pc_value)),
        poverview, width = overview_width, height = overview_height, dpi = 180, limitsize = FALSE)
      ggsave(file.path(plot_dir, sprintf("heritability_liability_overview_unfaceted_pc%d.pdf", pc_value)),
        poverview, width = overview_width, height = overview_height, limitsize = FALSE)
    }
    pdf(file.path(plot_dir, "heritability_liability_overview_unfaceted_all.pdf"),
      width = overview_width, height = overview_height, onefile = TRUE)
    for (poverview in overview_plots) print(poverview)
    dev.off()
  }
}

message("[DONE] Full tasks: ", nrow(totals), " | focus rows: ", nrow(focus), " | error/warning/not-run: ", nrow(errors))
message("[DONE] Summary directory: ", summary_dir)
