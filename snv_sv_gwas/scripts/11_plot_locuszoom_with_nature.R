#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

base <- "/path/to/EOSCZ_PROJECT"
ld_root <- file.path(base, "figure_analysis", "SV_SNV_LD", "LD_decay_public")
out_root <- file.path(ld_root, "locuszoom")
metadata_file <- file.path(out_root, "locuszoom_jobs.metadata.tsv")
nature_gwas <- Sys.getenv("EOSCZ_EXTERNAL_GWAS", unset = NA_character_)
window_bp <- 500000L
min_p <- 1e-300
y_cap <- 80

check_file <- function(path, label = "file") {
  if (length(path) != 1L || is.na(path) || !nzchar(as.character(path))) {
    stop(label, " path must be one non-missing character string; length=", length(path))
  }
  path <- as.character(path)
  if (!file.exists(path) || dir.exists(path) || file.info(path)$size == 0) {
    stop("Missing or empty ", label, ": ", path)
  }
}

norm_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x[toupper(x) == "X"] <- "23"
  x[toupper(x) == "Y"] <- "24"
  x
}

detect_col <- function(header, candidates, required = TRUE, label = "column") {
  idx <- match(tolower(candidates), tolower(header))
  idx <- idx[!is.na(idx)]
  if (length(idx)) return(header[idx[1L]])
  if (required) stop("Cannot find ", label, ". Available: ", paste(header, collapse = ", "))
  NA_character_
}

read_assoc <- function(path, label, require_id = TRUE) {
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  chr_col <- detect_col(header, c("CHR", "#CHR", "CHROM", "#CHROM", "chrom"), TRUE, paste(label, "CHR"))
  id_col <- detect_col(
    header,
    c("SNP", "ID", "MarkerName", "MARKER", "rsid", "RSID", "variant", "VARIANT", "varid"),
    require_id, paste(label, "variant ID")
  )
  pos_col <- detect_col(header, c("POS_T2T", "POS", "BP", "POSITION", "START"), TRUE, paste(label, "position"))
  p_col <- detect_col(header, c("P", "PVAL", "P_VALUE", "p", "pval"), TRUE, paste(label, "P"))
  select_cols <- unique(c(chr_col, id_col, pos_col, p_col))
  select_cols <- select_cols[!is.na(select_cols)]
  x <- fread(path, select = select_cols, showProgress = FALSE)
  setnames(x, chr_col, "CHR")
  setnames(x, pos_col, "POS")
  setnames(x, p_col, "P")
  if (!is.na(id_col)) setnames(x, id_col, "SNP")
  if (is.na(id_col)) x[, SNP := paste(norm_chr(CHR), POS, sep = ":")]
  x[, CHR := norm_chr(CHR)]
  x[, SNP := as.character(SNP)]
  x[, POS := suppressWarnings(as.integer(POS))]
  x[, P := suppressWarnings(as.numeric(P))]
  x <- x[!is.na(CHR) & !is.na(POS) & is.finite(P) & P > 0 & P <= 1]
  x[, LOGP := pmin(-log10(pmax(P, min_p)), y_cap)]
  x[, association_source := label]
  setkey(x, CHR, POS)
  x[]
}

assoc_cache <- new.env(parent = emptyenv())
get_assoc <- function(path, label, require_id = TRUE) {
  key <- paste(path, label, require_id, sep = "||")
  if (!exists(key, assoc_cache, inherits = FALSE)) {
    assign(key, read_assoc(path, label, require_id), assoc_cache)
  }
  copy(get(key, assoc_cache, inherits = FALSE))
}

read_ld <- function(path, lead_id) {
  check_file(path, "PLINK LD")
  x <- fread(path, showProgress = FALSE)
  required <- c("SNP_A", "SNP_B", "R2")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("LD file missing columns: ", paste(missing, collapse = ", "), " | ", path)
  x[, SNP_A := as.character(SNP_A)]
  x[, SNP_B := as.character(SNP_B)]
  x[, R2 := suppressWarnings(as.numeric(R2))]
  x <- x[is.finite(R2) & R2 >= 0 & R2 <= 1]
  map <- x[, .(R2 = max(R2)), by = .(SNP = SNP_B)]
  if (!lead_id %in% map$SNP) map <- rbind(map, data.table(SNP = lead_id, R2 = 1))
  map
}

ld_group <- function(r2) {
  out <- rep("LD unavailable", length(r2))
  out[!is.na(r2) & r2 < 0.2] <- "r2 < 0.2"
  out[!is.na(r2) & r2 >= 0.2 & r2 < 0.4] <- "r2 0.2-0.4"
  out[!is.na(r2) & r2 >= 0.4 & r2 < 0.6] <- "r2 0.4-0.6"
  out[!is.na(r2) & r2 >= 0.6 & r2 < 0.8] <- "r2 0.6-0.8"
  out[!is.na(r2) & r2 >= 0.8] <- "r2 >= 0.8"
  factor(out, levels = c(
    "LD unavailable", "r2 < 0.2", "r2 0.2-0.4", "r2 0.4-0.6",
    "r2 0.6-0.8", "r2 >= 0.8"
  ))
}

ld_colors <- c(
  "LD unavailable" = "#D9D9D9",
  "r2 < 0.2" = "#A6A6A6",
  "r2 0.2-0.4" = "#377EB8",
  "r2 0.4-0.6" = "#4DAF4A",
  "r2 0.6-0.8" = "#FF9F1C",
  "r2 >= 0.8" = "#E41A1C"
)

draw_one <- function(job, nature_dt) {
  chr <- norm_chr(job$lead_chr)
  pos <- as.integer(job$lead_pos)
  start <- max(1L, pos - window_bp)
  end <- pos + window_bp

  own <- get_assoc(job$regional_gwas_file, if (job$source_set %in% c("deepvariant", "paragraph")) "DeepVariant SGV" else job$source_set)
  own <- own[CHR == chr & POS >= start & POS <= end]
  if (!is.na(job$overlay_gwas_file) && nzchar(job$overlay_gwas_file)) {
    overlay <- get_assoc(job$overlay_gwas_file, "Paragraph SV")
    own <- rbindlist(list(own, overlay[CHR == chr & POS >= start & POS <= end]), fill = TRUE)
  }
  if (!nrow(own)) stop("No own-GWAS variants in window")
  if (!job$lead_id %in% own$SNP) stop("Lead ID is absent from regional GWAS points: ", job$lead_id)

  ld <- read_ld(job$ld_file, job$lead_id)
  own <- merge(own, ld, by = "SNP", all.x = TRUE, sort = FALSE)
  own[SNP == job$lead_id, R2 := 1]
  own[, LD_group := ld_group(R2)]
  own[, is_lead := SNP == job$lead_id]
  setorder(own, is_lead)
  nature <- nature_dt[CHR == chr & POS >= start & POS <= end]
  ymax <- max(c(own$LOGP, nature$LOGP), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  ymax <- min(y_cap, ymax * 1.08)
  xlim <- c(start, end) / 1e6

  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    layout(matrix(1:2, nrow = 2), heights = c(1, 1))
    par(mar = c(2.3, 4.4, 3.1, 1.0), mgp = c(2.4, 0.7, 0), tcl = -0.25)
    plot(NA, xlim = xlim, ylim = c(0, ymax), xlab = "", ylab = expression(-log[10](italic(P))),
         main = sprintf("Our GWAS | %s lead | %s | chr%s:%s", job$plot_class, job$source_set, chr, format(pos, big.mark = ",")),
         bty = "l", las = 1, cex.main = 0.9)
    abline(h = -log10(5e-8), lty = 2, col = "grey45")
    abline(h = -log10(1e-5), lty = 3, col = "grey65")
    for (grp in levels(own$LD_group)) {
      z <- own[LD_group == grp & is_lead == FALSE]
      if (nrow(z)) points(z$POS / 1e6, z$LOGP, pch = 20, cex = 0.52, col = ld_colors[[grp]])
    }
    lead <- own[is_lead == TRUE]
    points(lead$POS / 1e6, lead$LOGP, pch = 23, cex = 1.25, bg = "#7B3294", col = "black")
    legend("topright", legend = c("lead", names(ld_colors)), pch = c(23, rep(20, length(ld_colors))),
           pt.bg = c("#7B3294", rep(NA, length(ld_colors))), col = c("black", unname(ld_colors)),
           bty = "n", cex = 0.6, ncol = 2)

    par(mar = c(3.8, 4.4, 2.8, 1.0), mgp = c(2.4, 0.7, 0), tcl = -0.25)
    plot(NA, xlim = xlim, ylim = c(0, ymax), xlab = paste0("chr", chr, " position (Mb, T2T)"),
         ylab = expression(-log[10](italic(P))), main = "Nature primary_t2t | same locus",
         bty = "l", las = 1, cex.main = 0.9)
    abline(h = -log10(5e-8), lty = 2, col = "grey45")
    abline(h = -log10(1e-5), lty = 3, col = "grey65")
    abline(v = pos / 1e6, lty = 2, col = "#7B3294")
    if (nrow(nature)) points(nature$POS / 1e6, nature$LOGP, pch = 20, cex = 0.48, col = "grey30")
    if (!nrow(nature)) text(mean(xlim), ymax / 2, "No Nature variants in this window", col = "grey40")
  }

  pdf(paste0(job$plot_prefix, ".pdf"), width = 9, height = 7, useDingbats = FALSE)
  tryCatch(draw(), finally = dev.off())
  png(paste0(job$plot_prefix, ".png"), width = 2700, height = 2100, res = 300)
  tryCatch(draw(), finally = dev.off())

  own[, `:=`(
    lead_id = job$lead_id, lead_class = job$plot_class,
    lead_source_set = job$source_set, window_start = start, window_end = end
  )]
  fwrite(own, job$regional_table, sep = "\t", quote = FALSE, na = "NA")

  data.table(
    status = "PASS", plot_class = job$plot_class, source_set = job$source_set,
    lead_id = job$lead_id, lead_chr = chr, lead_pos = pos,
    n_own = nrow(own), n_ld_available = sum(!is.na(own$R2)),
    n_nature = nrow(nature), min_own_p = min(own$P),
    min_nature_p = if (nrow(nature)) min(nature$P) else NA_real_,
    regional_table = job$regional_table,
    plot_pdf = paste0(job$plot_prefix, ".pdf"),
    plot_png = paste0(job$plot_prefix, ".png"), error = NA_character_
  )
}

check_file(metadata_file, "locuszoom metadata")
check_file(nature_gwas, "Nature GWAS")
graphics.off()
jobs <- fread(metadata_file, showProgress = FALSE)
required_job_cols <- c(
  "plot_class", "source_set", "lead_id", "lead_chr", "lead_pos",
  "regional_gwas_file", "overlay_gwas_file", "ld_file",
  "plot_prefix", "regional_table"
)
missing_job_cols <- setdiff(required_job_cols, names(jobs))
if (length(missing_job_cols)) {
  stop(
    "Locuszoom metadata is missing required columns: ",
    paste(missing_job_cols, collapse = ", "),
    ". Rerun Step 13 with the current 10_prepare_locuszoom_ld.R."
  )
}
nature <- get_assoc(nature_gwas, "Nature primary_t2t", require_id = FALSE)

summary_rows <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i]
  message("[PLOT ", i, "/", nrow(jobs), "] ", job$plot_class, " | ", job$source_set, " | ", job$lead_id)
  summary_rows[[i]] <- tryCatch(
    draw_one(job, nature),
    error = function(e) data.table(
      status = "ERROR", plot_class = job$plot_class, source_set = job$source_set,
      lead_id = job$lead_id, lead_chr = job$lead_chr, lead_pos = job$lead_pos,
      error = conditionMessage(e)
    )
  )
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(out_root, "locuszoom_plot_summary.tsv"), sep = "\t", na = "NA")
for (cls in c("SGV", "SV")) {
  fwrite(summary_dt[plot_class == cls], file.path(out_root, cls, "locuszoom_plot_summary.tsv"), sep = "\t", na = "NA")
}
if (summary_dt[status == "ERROR", .N]) {
  stop(summary_dt[status == "ERROR", .N], " locuszoom plot(s) failed; see locuszoom_plot_summary.tsv")
}
message("[DONE] SGV and SV locuszoom plots and regional LD tables: ", out_root)
