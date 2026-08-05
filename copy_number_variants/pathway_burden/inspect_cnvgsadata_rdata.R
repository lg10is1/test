args_cli <- commandArgs(trailingOnly = TRUE)
get_cli_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args_cli[startsWith(args_cli, prefix)]
  if (length(hit) > 0) return(sub(prefix, "", hit[1], fixed = TRUE))
  default
}
default_project_root <- normalizePath(
  Sys.getenv("EOSCZ_PROJECT_ROOT", getwd()),
  winslash = "/",
  mustWork = FALSE
)
base <- get_cli_arg("burden-root", file.path(default_project_root, "cnv_analysis/pathway_burden"))
rdata <- file.path(base, "gene_set_reconstruction/cnvGSAdata_src/cnvGSAdata/data/gs_data_example.RData")
objs <- load(rdata)
cat("objects:", paste(objs, collapse=", "), "\n")
for (obj in objs) {
  x <- get(obj)
  cat("OBJECT", obj, "class", paste(class(x), collapse=","), "\n")
  str(x, max.level=2)
}

