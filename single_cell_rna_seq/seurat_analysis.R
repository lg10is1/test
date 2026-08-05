# PUBLIC-RELEASE METADATA-DRIVEN VERSION
# The source's embedded sample-to-phenotype mapping was replaced with a required de-identified TSV.
# Required columns: sample_id, time, mutation_group, disease_group, trio.
library(scCustomize)
library(future)
library(Seurat)
library(FLASHMM)
library(ggrastr)
library(patchwork)
library(data.table)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(gtable)
library(grid)



process_sample <- function(sample_dir) {
  data <- Read10X(data.dir = file.path(sample_dir, "filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(sample_dir), min.cells = 3, min.features = 200)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPL|^RPS")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
  seurat_obj <- NormalizeData(seurat_obj)
  #seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)  
  return(seurat_obj)
}

input_root <- Sys.getenv("EOSCZ_SCRNA_INPUT_DIR", unset = NA_character_)
metadata_path <- Sys.getenv("EOSCZ_SCRNA_SAMPLE_METADATA", unset = NA_character_)
output_dir <- Sys.getenv("EOSCZ_SCRNA_OUTPUT_DIR", unset = "seurat_output")
if (is.na(input_root) || !dir.exists(input_root)) stop("Set EOSCZ_SCRNA_INPUT_DIR to an authorized input directory")
if (is.na(metadata_path) || !file.exists(metadata_path)) stop("Set EOSCZ_SCRNA_SAMPLE_METADATA to a de-identified metadata TSV")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sample_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
sample_list <- lapply(sample_dirs, process_sample)
plan("multisession", workers = 15)
integrated_seurat <- Reduce(function(x, y) merge(x, y), sample_list)

#integrated_seurat[["RNA"]] <- JoinLayers(integrated_seurat[["RNA"]])




s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
s_genes <- intersect(s_genes, rownames(integrated_seurat))
g2m_genes <- intersect(g2m_genes, rownames(integrated_seurat))
integrated_seurat <- CellCycleScoring(
  integrated_seurat,
  s.features = s_genes,
  g2m.features = g2m_genes
)
#integrated_seurat[["RNA"]] <- split(integrated_seurat[["RNA"]], f = integrated_seurat$orig.ident)
integrated_seurat <- FindVariableFeatures(integrated_seurat, selection.method = "vst", nfeatures = 2000)
integrated_seurat <- ScaleData(integrated_seurat,vars.to.regress = c("percent.mt", "percent.ribo"))
integrated_seurat=RunPCA(integrated_seurat)
integrated_seurat <- IntegrateLayers(integrated_seurat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
save(integrated_seurat,file="integrated_seurat.RData")

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30, reduction = "harmony")
integrated_seurat <- FindClusters(integrated_seurat, resolution =0.4)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30, reduction = "harmony")
integrated_seurat[["RNA"]] <- JoinLayers(integrated_seurat[["RNA"]])


markers <- FindMarkers(integrated_seurat, ident.1 = "Exc", verbose = FALSE)
head(markers[which(markers[,2]>0),],n=20)


integrated_seurat=RenameIdents(integrated_seurat,c("0"="neuron","1"="deepExc","2"="dorsalNPC","3"="choroid_plexus","4"="mesenNPC","5"="rhombExc","7"="Cajal","8"="dorsalNPC_S",
"9"="IPC","10"="dorsalNPC_G2M","11"="cerebellar_NPC","13"="choroid_plexus","14"="dienNPC","15"="mesenNPC_S","16"="pericyte"))

cellvec=c("neuron","deepExc","rhombExc","dorsalNPC","ventralNPC","mesenNPC","dienNPC","cerebellarNPC","IPC","Cajal","choroidplexus")
marker=c("STMN2","DCX","BCL11B","TFAP2A","SOX2","VIM","EMX1","DLX2","OTX2","TCF7L2","UNC5C","EOMES","RELN","TTR")

exp=integrated_seurat@assays$RNA$data[marker,]
exp=as(exp,"matrix")
exp=t(exp)
exp=data.frame(exp)
exp$group=Idents(integrated_seurat)[rownames(exp)]
dhm=exp %>% group_by(group) %>% summarise(across(everything(), list(mean)))
dhm=data.frame(dhm)
colnames(dhm)=gsub("_1","",colnames(dhm))
rownames(dhm)=dhm[,1]
dhm=dhm[,-1]
dhm[]=apply(dhm,2,function(x){x/max(x)})
col_funb= colorRamp2(c(0,0.2,0.4,0.6,0.8,1),
c("white",brewer.pal(5,"Reds")))
dhm=dhm[cellvec,marker]
rownames(dhm)=c("neuron(level1)","deepExc(level4)","rhombExc(level4)","dorsalNPC(level3)","ventralNPC(level3)","mesenNPC(level2)","dienNPC(level2)","cerebellarNPC(level3)",
"IPC","Cajal","choroidplexus(level1)")
f1d=Heatmap(as.matrix(dhm),show_row_dend = FALSE,name="exp",
show_column_dend=FALSE,col=col_funb,
row_title_rot = 0,row_order=rownames(dhm),column_order=marker,show_heatmap_legend=F)

pdf("markerhm.pdf",height=5,width=6)
draw(f1d)
dev.off()


sample_metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
required_metadata <- c("sample_id", "time", "mutation_group", "disease_group", "trio")
missing_metadata <- setdiff(required_metadata, colnames(sample_metadata))
if (length(missing_metadata) > 0) stop(paste("Missing metadata columns:", paste(missing_metadata, collapse = ", ")))
metadata_index <- match(integrated_seurat@meta.data$orig.ident, sample_metadata$sample_id)
if (anyNA(metadata_index)) stop("Every orig.ident value must have one de-identified metadata row")
integrated_seurat@meta.data$time <- sample_metadata$time[metadata_index]
integrated_seurat@meta.data$mutation <- sample_metadata$mutation_group[metadata_index]
integrated_seurat@meta.data$SCZ <- sample_metadata$disease_group[metadata_index]
integrated_seurat@meta.data$trio <- sample_metadata$trio[metadata_index]
table(integrated_seurat@meta.data[,c("orig.ident","cluster")])
integrated_seurat@meta.data$cluster=Idents(integrated_seurat)


integrated_seurat=subset(integrated_seurat,new_cluster != "pericyte")
p_time=DimPlot(integrated_seurat, reduction = "umap",group.by = "time",raster=T)+theme(
  axis.text = element_blank(),                    
  axis.ticks = element_blank(),                 
  axis.title = element_text(vjust = 0,size = 8),             
  legend.key.size = unit(0.5, "lines"),  
  legend.text = element_text(size = 12),  
  legend.spacing = unit(0.1, "cm"),      
  legend.margin = margin(0, 0, 0, 0),    
  legend.box.spacing = unit(0.1, "cm"),  
  legend.key.height = unit(0.3, "cm"),  
  legend.key.width = unit(0.3, "cm")
)+scale_color_brewer(palette = "Set1")



polychrome_colors <- DiscretePalette_scCustomize(11, palette = "polychrome")

p_cluster=DimPlot(integrated_seurat, reduction = "umap",group.by = "level1", cols = polychrome_colors,raster=T)+theme(
  axis.text = element_blank(),                    
  axis.ticks = element_blank(),                 
  axis.title = element_blank(),             
  legend.key.size = unit(0.5, "lines"),  
  legend.text = element_text(size = 12),  
  legend.spacing = unit(0.1, "cm"),      
  legend.margin = margin(0, 0, 0, 0),    
  legend.box.spacing = unit(0.1, "cm"),  
  legend.key.height = unit(0.3, "cm"),  
  legend.key.width = unit(0.3, "cm") 
)

point_data <- data.frame(
  x = 1:13,
  y = rep(-0.1, 13),
  color=factor(1:13)
)

pdf("plot_up.pdf",height=4,width=10)
p_time+p_cluster+plot_layout(guides="collect")
dev.off()

pdf("level1_umap_tel.pdf",height=5,width=5)
DimPlot(integrated_seurat, reduction = "umap",group.by = "level1", cols = polychrome_colors,raster=T)+theme(
  axis.text = element_blank(),                    
  axis.ticks = element_blank(),                 
  axis.title = element_blank(),             
  legend.key.size = unit(0.5, "lines"),  
  legend.text = element_text(size = 12),  
  legend.spacing = unit(0.1, "cm"),      
  legend.margin = margin(0, 0, 0, 0),    
  legend.box.spacing = unit(0.1, "cm"),  
  legend.key.height = unit(0.3, "cm"),  
  legend.key.width = unit(0.3, "cm") 
)
dev.off()

meta=integrated_seurat@meta.data
meta=meta[which(meta$new_cluster!="pericyte"),]

prop=as.matrix(table(meta$time,meta$new_cluster))
prop=t(prop)
prop[]=apply(prop,2,function(x){x/sum(x)})
prop=prop[order(rowSums(prop),decreasing=T),]
col_fun <- colorRamp2(c(0, 0.4), c("white", "darkred"))
pdf("prop.pdf",height=5,width=5)
Heatmap(prop, name = "proportion",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        column_names_side = "bottom",
)
dev.off()



expression_matrix <- GetAssayData(integrated_seurat, layer = "data") 
d=expression_matrix[c("SLC17A7","GAD1","GFAP","TBR1","EOMES","SATB2","RELN","CDH12","LHFPL4","GLIS3"),]
d=t(as.matrix(d))
d=data.frame(d)
d$cell=rownames(d)
long_df <- data.frame(d) %>%
  tidyr::pivot_longer(
    cols = "SLC17A7":"GLIS3",
    names_to = "gene",
    values_to = "exp"
)
umap_coords <- Embeddings(integrated_seurat, reduction = "umap")
umap1 <- umap_coords[, 1]  
umap2 <- umap_coords[, 2]

draw=data.frame(data.frame(long_df),umap_coords[match(long_df$cell,rownames(umap_coords)),1:2],time=meta[match(long_df$cell,rownames(meta)),"time"])

p15=ggplot(draw, aes(x = umap1, y = umap2, color = SLC17A6)) +
  geom_point(size = 0.01)+facet_grid(time~sample)+
  scale_color_continuous(limits = c(0, 3),low="grey", high="#272973") +theme_classic()+
  theme(legend.position = "none",axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title=element_text(color = "black", size = 6),
            strip.background = element_rect(fill = NA, color = NA))
            
            
for (trioidx in 1:3){
int2 <- subset(integrated_seurat,trio == trioidx)
for (i in unique(integrated_seurat$new_cluster)) {
  
  int1 <- subset(int2, new_cluster == i)
  
  for (t in c(30, 60,120)) {
    
    if (sum(int1$time == t) < 5) next
    
    int <- subset(int1, time == t)
    
    if (min(table(int$SCZ)) < 5 | length(table(int$SCZ)) < 2) next
    

    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
    counts_mat <- GetAssayData(int, assay = "RNA", layer = "counts")
    counts_mat <- counts_mat[rowSums(counts_mat) >= 50, ]
    Y <- log2(1 + counts_mat)
    meta <- int@meta.data
    
# Source comment removed during English public-release translation.
    meta$SCZ <- factor(meta$SCZ, levels = c("control", "case"))  
    meta$orig.ident <- as.factor(meta$orig.ident)                
    meta$libsize <- colSums(counts_mat)
    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
# Source comment removed during English public-release translation.
    X <- model.matrix(~ log(libsize) + SCZ, data = meta)
# Source comment removed during English public-release translation.
    Z <- model.matrix(~ orig.ident, data = meta)
    d <- ncol(Z)
    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
    tryCatch({
      Y_mat <- as.matrix(Y)
      n <- nrow(X)
      
      XX    <- crossprod(X)
      XY    <- crossprod(X, t(Y_mat))
      ZX    <- crossprod(Z, X)
      ZY    <- crossprod(Z, t(Y_mat))
      ZZ    <- crossprod(Z)
      Ynorm <- rowSums(Y_mat * Y_mat)
      
      rm(Y_mat, Y, counts_mat); gc()
      
# Source comment removed during English public-release translation.
      fit <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)
      
      # ============================================================
# Source comment removed during English public-release translation.
      # ============================================================
      ct <- numeric(ncol(X))
      idx_scz <- which(colnames(X) == "SCZcase")  
      ct[idx_scz] <- 1
      
      test_res <- as.data.frame(lmmtest(fit, contrast = ct))
      
      if (nrow(test_res) > 0) {
        markers <- data.frame(
          gene     = rownames(test_res),
          avg_log2FC = test_res$`_coef`,
          p_val    = test_res$`_p`,
          p_val_adj = p.adjust(test_res$`_p`, method = "BH"),
          cluster  = i,
          time     = t,
          trio=trioidx,
          stringsAsFactors = FALSE
        )
        fwrite(markers,"lmm_trio.csv",append=T)
        print(c(i, t, trioidx, nrow(markers)))
      }
      
      rm(fit, XX, XY, ZX, ZY, ZZ, Ynorm); gc()
      
    }, error = function(e) {
      print(c(i, t, paste("ERROR:", e$message)))
    })
  }
}
}


for (i in unique(integrated_seurat$new_cluster)) {
  
  int1 <- subset(integrated_seurat, new_cluster == i)
  
  for (t in c(15,30,60,120)) {
    
    if (sum(int1$time == t) < 5) next
    
    int <- subset(int1, time == t)
    
    if (min(table(int$SCZ)) < 5 | length(table(int$SCZ)) < 2) next
    

    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
    counts_mat <- GetAssayData(int, assay = "RNA", layer = "counts")
    counts_mat <- counts_mat[rowSums(counts_mat) >= 50, ]
    Y <- log2(1 + counts_mat)
    meta <- int@meta.data
    
# Source comment removed during English public-release translation.
    meta$SCZ <- factor(meta$SCZ, levels = c("control", "case"))  
    meta$orig.ident <- as.factor(meta$orig.ident)                
    meta$libsize <- colSums(counts_mat)
    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
# Source comment removed during English public-release translation.
    X <- model.matrix(~ log(libsize) + SCZ, data = meta)
# Source comment removed during English public-release translation.
    Z <- model.matrix(~ orig.ident, data = meta)
    d <- ncol(Z)
    
    # ============================================================
# Source comment removed during English public-release translation.
    # ============================================================
    tryCatch({
      Y_mat <- as.matrix(Y)
      n <- nrow(X)
      
      XX    <- crossprod(X)
      XY    <- crossprod(X, t(Y_mat))
      ZX    <- crossprod(Z, X)
      ZY    <- crossprod(Z, t(Y_mat))
      ZZ    <- crossprod(Z)
      Ynorm <- rowSums(Y_mat * Y_mat)
      
      rm(Y_mat, Y, counts_mat); gc()
      
# Source comment removed during English public-release translation.
      fit <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d)
      
      # ============================================================
# Source comment removed during English public-release translation.
      # ============================================================
      ct <- numeric(ncol(X))
      idx_scz <- which(colnames(X) == "SCZcase")  
      ct[idx_scz] <- 1
      
      test_res <- as.data.frame(lmmtest(fit, contrast = ct))
      
      if (nrow(test_res) > 0) {
        markers <- data.frame(
          gene     = rownames(test_res),
          avg_log2FC = test_res$`_coef`,
          p_val    = test_res$`_p`,
          p_val_adj = p.adjust(test_res$`_p`, method = "BH"),
          cluster  = i,
          time     = t,
          stringsAsFactors = FALSE
        )
        fwrite(markers,"lmm.csv",append=T)
        print(c(i, t,nrow(markers)))
      }
      
      rm(fit, XX, XY, ZX, ZY, ZZ, Ynorm); gc()
      
    }, error = function(e) {
      print(c(i, t, paste("ERROR:", e$message)))
    })
  }
}
