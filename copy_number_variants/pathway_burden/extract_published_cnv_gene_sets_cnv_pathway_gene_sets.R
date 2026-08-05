# Dependencies must be installed separately; this script performs no package installation.
#[0]----
# ## Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
# 
#   "clusterProfiler",
#   "org.Hs.eg.db",
#   "KEGGREST",
#   "ReactomePA",
#   "biomaRt"
# ))
# 
# options(timeout = 1000)
# ## CRAN
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(ReactomePA)
library(biomaRt)
library(dplyr)
library(purrr)

#[1]----GO-derived gene sets (Neurof_Go*)
get_go_genes <- function(go_ids) {
  genes <- map(go_ids, function(go) {
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go,
        keytype = "GO",
        columns = "SYMBOL"
      )$SYMBOL
    )
  }) %>% unlist() %>% unique() %>% na.omit()
  return(genes)
}

Neurof_GoNeuronBody <- get_go_genes("GO:0043025")
length(Neurof_GoNeuronBody)  # ~309

Neurof_GoSynaptic <- get_go_genes(c(
  "GO:0045202",  # synapse
  "GO:0050808"   # synapse organization
))

Neurof_GoNeuronProj <- get_go_genes(c(
  "GO:0043005",  # neuron projection
  "GO:0031175"   # neuron projection development
))

Neurof_GoNervTransm <- get_go_genes(c(
  "GO:0019226",  # transmission of nerve impulse
  "GO:0007268"   # synaptic transmission
))

Neurof_GoNervSysDev <- get_go_genes("GO:0007399")
Neurof_GoNervSysDev_CNS <- get_go_genes("GO:0007417")

#[2]----KEGG synaptic pathways
kegg_ids <- c(
  "hsa04725", "hsa04724", "hsa04728",
  "hsa04727", "hsa04726", "hsa04721",
  "hsa04723", "hsa04720", "hsa04730"
)

get_kegg_genes <- function(pid) {
  k <- keggGet(pid)[[1]]$GENE
  if (is.null(k)) return(NULL)
  k[seq(2, length(k), 2)] %>%
    gsub(";.*", "", .)
}

Neurof_KeggSynaptic <- map(kegg_ids, get_kegg_genes) %>%
  unlist() %>% unique()

#[3]----Axon guidance pathways (Neurof_PathwaysAxonG; approximate reconstruction)
## KEGG Axon guidance
axon_kegg <- get_kegg_genes("hsa04360")

## Reactome Axon-related
axon_reactome <- enrichPathway(
  gene = keys(org.Hs.eg.db, "ENTREZID"),
  organism = "human",
  readable = TRUE
)@result %>%
  filter(grepl("axon|NCAM|netrin|reelin", Description, ignore.case = TRUE)) %>%
  pull(geneID) %>%
  paste(collapse = "/") %>%
  strsplit("/") %>%
  unlist()

Neurof_PathwaysAxonG <- unique(c(axon_kegg, axon_reactome))

#[4]----Union gene sets(Neurof_Union*)
Neurof_UnionInclusive <- unique(c(
  Neurof_KeggSynaptic,
  Neurof_GoNervTransm,
  Neurof_GoNeuronProj,
  Neurof_GoNeuronBody,
  Neurof_GoSynaptic,
  Neurof_GoNervSysDev,
  Neurof_PathwaysAxonG
))

#[5]----Stringent union gene sets (appearing in two or more sources)
all_sets <- list(
  Neurof_KeggSynaptic,
  Neurof_GoNervTransm,
  Neurof_GoNeuronProj,
  Neurof_GoNeuronBody,
  Neurof_GoSynaptic,
  Neurof_GoNervSysDev,
  Neurof_PathwaysAxonG
)

Neurof_UnionStringent <- unlist(all_sets) %>%
  table() %>%
  .[. >= 2] %>%
  names()

#[6]----

