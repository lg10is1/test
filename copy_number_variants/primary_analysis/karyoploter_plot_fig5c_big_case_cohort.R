suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- normalizePath(sub(file_arg, "", args[grep(file_arg, args)]), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_path)

input_table <- file.path(script_dir, "gene_frequencies_with_chr_length_case_cohort.txt")
output_pdf <- file.path(script_dir, "Fig5Cbig_CHM13_karyotype_plot_case_cohort.pdf")
output_png <- file.path(script_dir, "Fig5Cbig_CHM13_karyotype_plot_case_cohort.png")
censat_bed <- file.path(script_dir, "chm13v2.0_censat_v2.1.bed")
sd_bed <- file.path(script_dir, "chm13v2.0_SD.bed")
telomere_bed <- file.path(script_dir, "chm13v2.0_telomere.bed")

autosomes <- paste0("chr", 1:22)
df <- read.delim(input_table, header = TRUE, sep = "\t")
df <- df[df$Chromosome %in% autosomes, ]

chm13_chrom_sizes <- toGRanges(data.frame(
  chr = autosomes,
  start = rep(1, length(autosomes)),
  end = c(
    248387328, 242696752, 201105948, 193574945, 182045439, 172126628,
    160567428, 146259331, 150617247, 134758134, 135127769, 133324548,
    113566686, 101161492, 99753195, 96330374, 84276897, 80542538,
    61707364, 66210255, 45090682, 51324926
  )
))

cen_data <- read.delim(censat_bed, header = FALSE, sep = "\t")
colnames(cen_data) <- c("Chromosome", "Start", "End", "Name")
cen_data <- cen_data[cen_data$Chromosome %in% autosomes, ]

sd_data <- read.delim(sd_bed, header = FALSE, sep = "\t")
colnames(sd_data) <- c("Chromosome", "Start", "End", "Name")
sd_data <- sd_data[sd_data$Chromosome %in% autosomes, ]

tel_data <- read.delim(telomere_bed, header = FALSE, sep = "\t")
colnames(tel_data) <- c("Chromosome", "Start", "End")
tel_data <- tel_data[tel_data$Chromosome %in% autosomes, ]

pdf(output_pdf, width = 10, height = 8)
kp <- plotKaryotype(genome = chm13_chrom_sizes, plot.type = 6, chromosomes = "all")
kpDataBackground(kp, data.panel = 1, col = "#EEEEEE")
kpRect(kp, chr = sd_data$Chromosome, x0 = sd_data$Start, x1 = sd_data$End,
       y0 = 0.03, y1 = 0.97, col = "lightblue", border = NA)
kpRect(kp, chr = cen_data$Chromosome, x0 = cen_data$Start, x1 = cen_data$End,
       y0 = 0.03, y1 = 0.97, col = "lightgrey", border = NA)
kpPoints(kp, chr = df$Chromosome, x = (df$Start + df$End) / 2, y = 0.5,
         col = "#da7271", bg = "#da7271", cex = 1, pch = 25)
dev.off()

png(output_png, width = 2000, height = 1600, res = 200)
kp <- plotKaryotype(genome = chm13_chrom_sizes, plot.type = 6, chromosomes = "all")
kpDataBackground(kp, data.panel = 1, col = "#EEEEEE")
kpRect(kp, chr = sd_data$Chromosome, x0 = sd_data$Start, x1 = sd_data$End,
       y0 = 0.03, y1 = 0.97, col = "lightblue", border = NA)
kpRect(kp, chr = cen_data$Chromosome, x0 = cen_data$Start, x1 = cen_data$End,
       y0 = 0.03, y1 = 0.97, col = "lightgrey", border = NA)
kpPoints(kp, chr = df$Chromosome, x = (df$Start + df$End) / 2, y = 0.5,
         col = "#da7271", bg = "#da7271", cex = 1, pch = 25)
dev.off()

cat("Saved:", output_pdf, "\n")
cat("Saved:", output_png, "\n")
