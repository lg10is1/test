library(data.table)
library(qqman)
library(ggsci)
library(ggplot2)
library(dplyr)
library(cowplot)
# ============================================================
# 0. Load GWAS results.
#    sgv_sites contains SNVs and indels shorter than 50 bp.
#    sv_sites contains structural variants of at least 50 bp.
# ============================================================

in_dir <- "/path/to/EOSCZ_PROJECT/figure_analysis/01.GWAS_figure.public"
deepvariant_paragraph_dir <- "/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph"

deepvariant_snv_gwas_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/deepvar/03_gwas/SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA"
)

paragraph_sv_gwas_file <- file.path(
  deepvariant_paragraph_dir,
  "deepvar_gwas/paragraph_test/03_gwas/SCZ.paragraph_test.pcsrc_deepvar_pc20_grm_deepvar_with_batch.mlm.geno0.1.maf0.01.fastGWA"
)

sgv_sites <- fread(
  file.path(
    in_dir,
    "all_sets.merged_minP.SNV_INDEL_lt50bp.remove_chr_edge_1Mb.standardized.tsv"
  )
)

sv_sites <- fread(
  file.path(
    in_dir,
    "all_sets.merged_minP.SV_ge50bp.remove_chr_edge_1Mb.standardized.tsv"
  )
)

# ============================================================
# 1. Validate required columns.
# ============================================================

required_cols <- c("CHR", "POS", "min_P")

missing_sgv <- setdiff(required_cols, names(sgv_sites))
missing_sv <- setdiff(required_cols, names(sv_sites))

if (length(missing_sgv) > 0) {
  stop("sgv_sites is missing required columns: ", paste(missing_sgv, collapse = ", "))
}

if (length(missing_sv) > 0) {
  stop("sv_sites is missing required columns: ", paste(missing_sv, collapse = ", "))
}

message("Loaded sgv_sites: ", nrow(sgv_sites))
message("Loaded sv_sites: ", nrow(sv_sites))



sgv_sites$P=sgv_sites$min_P
sv_sites$P=sv_sites$min_P

sgv_sites$CHR=as.numeric(sgv_sites$CHR)
sv_sites$CHR=as.numeric(sv_sites$CHR)

qq(sv_sites$P)

standardize_external_gwas_for_plot <- function(file, type, source) {
  if (!file.exists(file)) {
    stop("GWAS file does not exist: ", file)
  }
  
  x <- fread(file)
  
  if ("#CHROM" %in% names(x) && !"CHR" %in% names(x)) {
    setnames(x, "#CHROM", "CHR")
  }
  
  if ("CHROM" %in% names(x) && !"CHR" %in% names(x)) {
    setnames(x, "CHROM", "CHR")
  }
  
  if ("BP" %in% names(x) && !"POS" %in% names(x)) {
    x[, POS := BP]
  }
  
  required_cols_external <- c("CHR", "POS", "P")
  missing_external <- setdiff(required_cols_external, names(x))
  
  if (length(missing_external) > 0) {
    stop(
      "External GWAS file missing columns: ",
      paste(missing_external, collapse = ", "),
      "\nfile: ",
      file
    )
  }
  
  x[, CHR := as.numeric(sub("^chr", "", as.character(CHR), ignore.case = TRUE))]
  x[, POS := as.numeric(POS)]
  x[, P := as.numeric(P)]
  x <- x[!is.na(CHR) & !is.na(POS) & !is.na(P) & P > 0 & P <= 1]
  x[, type := type]
  x[, source := source]
  
  x[, .(CHR, POS, P, type, source)]
}





snv2=sgv_sites[, .(CHR, POS, P)]
snv2=snv2[order(snv2$CHR,snv2$POS),]
snv2$type='snv'
snv2$source='pangenie'
sv=sv_sites[, .(CHR, POS, P)]
sv$type='sv'
sv$source='pangenie'

deepvariant_snv <- standardize_external_gwas_for_plot(
  file = deepvariant_snv_gwas_file,
  type = "snv",
  source = "deepvariant"
)

paragraph_sv <- standardize_external_gwas_for_plot(
  file = paragraph_sv_gwas_file,
  type = "sv",
  source = "paragraph"
)

snv2=filter(snv2,CHR %in% 1:22)
sv=filter(sv,CHR %in% 1:22)
deepvariant_snv=filter(deepvariant_snv,CHR %in% 1:22)
paragraph_sv=filter(paragraph_sv,CHR %in% 1:22)

gwas_all=rbindlist(list(snv2,deepvariant_snv,sv,paragraph_sv),use.names=TRUE,fill=TRUE)
gwas_all=gwas_all[!is.na(P) & P > 0 & P <= 1]
gwas_all=gwas_all[order(gwas_all$CHR,gwas_all$POS),]
gwasResults =data.frame(gwas_all)
ymax <- ceiling(max(-log10(gwasResults$P), na.rm = TRUE) + 2)

chr_len <- gwasResults %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(POS))
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len)
Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, POS) %>%
  mutate( BPcum = POS + total)
X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) +min(BPcum) ) / 2 )

snv_plot=filter(Snp_pos,type=='snv')
sv_plot=filter(Snp_pos,type=='sv')

p_sv <- ggplot(sv_plot, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1, stroke=0) +
  #      
  scale_color_manual(values = rep(c("#47a1a2", "#272973"), 22 )) +
  scale_x_continuous(labels = X_axis$CHR, breaks = X_axis$center) +
  #         X      gap
  scale_y_continuous(expand = c(0, 0),limits = c(0,ymax) ) + 
  #        
  geom_hline(yintercept = c(-log10(5e-6), -log10(5e-8)), color = c('grey', 'red'), linetype = c("dotted", "twodash")) + 
  #      
  theme_cowplot(8) +theme(legend.position="none")+labs(x="CHR")



p_snv2 <- ggplot(snv_plot, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1, stroke=0) +
  #      
  scale_color_manual(values = rep(c("#47a1a2", "#272973"), 22 )) +
  scale_x_continuous(labels = X_axis$CHR, breaks = X_axis$center) +
  #         X      gap
  scale_y_continuous(expand = c(0, 0),limits = c(0,ymax) ) + 
  #        
  geom_hline(yintercept = c(-log10(5e-6), -log10(5e-8)), color = c('grey', 'red'), linetype = c("dotted", "twodash")) + 
  #      
  theme_cowplot(8) +theme(legend.position="none")+labs(x="CHR")

p_snv <- ggplot(snv_plot, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1, stroke=0) +
  #      
  scale_color_manual(values = rep(c("#47a1a2", "#272973"), 22 )) +
  scale_x_continuous(labels = X_axis$CHR, breaks = X_axis$center, position = "top") +
  #         X      gap
  scale_y_reverse(expand = c(0, 0),limits = c(ymax,0)) +
  #        
  geom_hline(yintercept = c(-log10(5e-6), -log10(5e-8)), color = c('grey', 'red'), linetype = c("dotted", "twodash")) + 
  #      
  theme_cowplot(8) +theme(legend.position="none",axis.text.x = element_blank(),  #    x      
                          axis.ticks.x = element_blank())+labs(x=NULL)

p2=plot_grid(p_sv,p_snv,nrow=2)

#ggsave(p_snv2,file="/path/to/EOSCZ_PROJECT/Figure3/snv_gwas.tiff",width=1889,height=1299,units='px',dpi=300,bg='white')
ggsave(p_sv,file="/path/to/EOSCZ_PROJECT/figure_analysis/sv_gwas.public.tiff",width=1889,height=1299,units='px',dpi=300,bg='white')
ggsave(p2,file="/path/to/EOSCZ_PROJECT/figure_analysis/sv_snv_gwas.public.tiff",width=155.1,height=106.7,units='mm',dpi=300,bg='white')

