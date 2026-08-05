# Load necessary libraries
x <- commandArgs(trailingOnly = TRUE)[1]
x=as.numeric(x)

library(data.table)
library(logistf)
#to enable parallel tasks, TR matrix was split into smaller chunks
#NGS GWAS and TGS GWAS used similar codes on different TR matrix
d3=fread(paste("tmp/dscz/",x,sep=""),data.table=F)
d0=fread(paste("tmp/dcontrol/",x,sep=""),data.table=F)


load("lord.RData") #TRID for each chunk
lord=lord[[x]]
info=fread("sample.info",data.table=F)
colnames(info)[2]="SCZ"

n_SCZ_haps=ncol(d3)
n_SCZ_sample=n_SCZ_haps/2
n_control_haps=ncol(d0)
n_control_sample=n_control_haps/2
for (i in 1:nrow(d0)) {
  ordid <- as.character(lord[i, 1])
  sl <- unlist(d3[i, ])
  ssd <- unlist(sczsd[i, ])
  cl <- unlist(d0[i, ])
  sl[which(ssd == 1)] <- NA
  cl <- (cl[1:n_control_sample] + cl[(n_control_sample+1):n_control_haps]) / 2
  SL <- (sl[1:n_SCZ_sample] + sl[(n_SCZ_sample+1):n_SCZ_haps]) / 2

  cmax <- max(na.omit(cl))
  cmin <- min(na.omit(cl))
  smax <- max(na.omit(SL))
  smin <- min(na.omit(SL))
  
  if (cmax == cmin | smax == smin) {next}
  
  tr <- c(SL, cl)
  tr <- qnorm((rank(tr, na.last = "keep") - 0.5) / sum(!is.na(tr)))
  tr <- data.frame(tr = tr, info[,-1, drop = FALSE])
  tr <- na.omit(tr)
  n <- nrow(tr)
  
  skip_flag <- FALSE
  
  model <- tryCatch(
    {
      glm(SCZ ~ ., data = tr, family = "binomial")
    },
    error = function(e) {
      skip_flag <<- TRUE
      return(NULL)
    }
  )
  
  if (skip_flag || is.null(model)) {
    next
  }
  
  int <- summary(model)$coefficients
  tr_row <- match("tr", rownames(int))
  
  if (!is.na(tr_row)) {
    rescon <- c(
      id = ordid, 
      unlist(int[tr_row, ]), 
      n = n, 
      cmax = cmax, 
      cmin = cmin, 
      smax = smax, 
      smin = smin
    )
    
    fwrite(as.list(rescon), "trlogistic_tgs", append = TRUE, sep = "\t")
  }
}


