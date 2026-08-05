required_packages <- c("dbscan", "parallel", "Matrix","matrixStats","data.table","tidyr","dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

#MC_* are TR count matrix (MC field of trgt output) for case and control. each row has $TRID_$motiforder as column 1 (using awk command). trgt output ensure same row number and order, allowing direct combination by columns.
d_scz=fread("MC_SCZ",data.table=F)
d_control=fread("MC_control",data.table=F)
rownames(d_scz)=d_scz[,1]
d_scz=d_scz[,-1]
rownames(d_control)=d_control[,1]
d_control=d_control[,-1]
sczmax=rowMaxs(d_scz)
sczmin=rowMins(d_scz)
controlmax=rowMaxs(d_control)
controlmin=rowMins(d_control)
d=cbind(d_control,d_scz)

detect_outliers_safe <- function(x, row_idx, minPts = 5, eps_quantile = 0.9, min_eps = 5.0) {
 
  n_samples <- length(x)
  
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1) {
    return(NULL)  
  }
  
  na_idx <- is.na(x)
  x_non_na <- x[!na_idx]
  non_na_pos <- which(!na_idx)  
  
  x_sorted <- sort(x_non_na)
  diffs <- diff(x_sorted)
  eps <- max(min_eps, quantile(diffs, probs = eps_quantile, na.rm = TRUE), na.rm = TRUE)
  
  x_mat <- matrix(x_non_na, ncol = 1)
  db <- dbscan(x_mat, eps = eps, minPts = minPts)
  
  if (all(db$cluster == 0)) {
    return(NULL)
  }
  
  cluster_counts <- table(db$cluster[db$cluster > 0])
  main_cluster <- as.integer(names(cluster_counts)[which.max(cluster_counts)])
  main_center <- median(x_non_na[db$cluster == main_cluster])
  
  noise_idx <- which(db$cluster == 0)
  if (length(noise_idx) == 0) {
    return(NULL)  
  }
  
  noise_vals <- x_non_na[noise_idx]
  higher <- noise_vals > main_center
  
  triplets <- data.frame(
    i = row_idx,
    j = non_na_pos[noise_idx],
    val = ifelse(higher, 1L, -1L)
  )
  
  return(triplets)
}



minPts <- 5
eps_quantile <- 0.9       
min_eps <- 10.0   
#n_cores <- parallel::detectCores() - 1  
n_cores=19
n_loci <- nrow(d)
n_samples <- ncol(d)
locus_ids <- rownames(d)
sample_ids <- colnames(d)

results_list <- mclapply(seq_len(n_loci), function(i) {
  if (i %% 1000 == 0) cat(sprintf("Analyze %d / %d \n", i, n_loci))
  
  row_vals <- d[i, ]
  detect_outliers_safe(row_vals, i, minPts, eps_quantile, min_eps)
}, mc.cores = n_cores, mc.preschedule = TRUE)

triplets <- do.call(rbind, Filter(Negate(is.null), results_list))
save(triplets,file="triplets.RData")
sparse_mat <- sparseMatrix(
  i = triplets$i,
  j = triplets$j,
  x = triplets$val,
  dims = c(n_loci, n_samples)  
)

locusid=rownames(d)
hapid=colnames(d)

spcontrol <- subset(triplets, triplets[, 2] <= ncol(d_control))
spcase <- subset(triplets, triplets[, 2] > ncol(d_control))

res=data.frame(ID=rownames(d_scz),ID1=rownames(d_scz)) 
res=separate(res,ID1,into=c("posid","motifid"),sep=":")
res$mintr_case=sczmin
res$maxtr_case=sczmax
res$mintr_control=controlmin
res$maxtr_control=controlmax

res$tre_control=0
res$tre_case=0
res$trc_control=0
res$trc_case=0

tre_case=table(spcase[which(spcase[,3]==1),1])
res[names(tre_case),"tre_case"]=tre_case

trc_case=table(spcase[which(spcase[,3]==(-1)),1])
res[names(trc_case),"trc_case"]=trc_case

tre_control=table(spcontrol[which(spcontrol[,3]==1),1])
res[names(tre_control),"tre_control"]=tre_control

trc_control=table(spcontrol[which(spcontrol[,3]==(-1)),1])
res[names(trc_control),"trc_control"]=trc_control

annot=fread("trannot.txt",data.table=F)
res$pos=annot$pos[match(res$posid,annot$posid)]
res$cre=annot$cre[match(res$posid,annot$posid)]

burden=data.frame()
for(cre in unique(res$cre)){
  if(is.na(cre)){next}
  int=res[which(res$cre==cre),]
  trecon=sum(int$tre_control)
  trecase=sum(int$tre_case)
  model=binom.test(trecase, trecase+trecon, p = sum(res$tre_case,na.rm=T)/(sum(res$tre_case,na.rm=T)+sum(res$tre_control,na.rm=T)), alternative = "greater")
  p=model$p.value
  RR=model$estimate/(sum(res$tre_case,na.rm=T)/(sum(res$tre_case,na.rm=T)+sum(res$tre_control,na.rm=T)))
  burden=rbind(burden,data.frame(cre=cre,type="expansion",p=p,RR=RR))
  
  trccon=sum(int$trc_control)
  trccase=sum(int$trc_case)
  model=binom.test(trccase, trccase+trccon, p = sum(res$trc_case,na.rm=T)/(sum(res$trc_case,na.rm=T)+sum(res$trc_control,na.rm=T)), alternative = "greater")
  p=model$p.value
  RR=model$estimate/(sum(res$trc_case,na.rm=T)/(sum(res$trc_case,na.rm=T)+sum(res$trc_control,na.rm=T)))
  burden=rbind(burden,data.frame(cre=cre,type="contraction",p=p,RR=RR))  
}

for(cre in unique(res$pos)){
  if(is.na(cre)){next}
  int=res[which(res$pos==cre),]
  trecon=sum(int$tre_control)
  trecase=sum(int$tre_case)
  model=binom.test(trecase, trecase+trecon, p = sum(res$tre_case,na.rm=T)/(sum(res$tre_case,na.rm=T)+sum(res$tre_control,na.rm=T)), alternative = "greater")
  p=model$p.value
  RR=model$estimate/(sum(res$tre_case,na.rm=T)/(sum(res$tre_case,na.rm=T)+sum(res$tre_control,na.rm=T)))
  burden=rbind(burden,data.frame(cre=cre,type="expansion",p=p,RR=RR))
  
  trccon=sum(int$trc_control)
  trccase=sum(int$trc_case)
  model=binom.test(trccase, trccase+trccon, p = sum(res$trc_case,na.rm=T)/(sum(res$trc_case,na.rm=T)+sum(res$trc_control,na.rm=T)), alternative = "greater")
  p=model$p.value
  RR=model$estimate/(sum(res$trc_case,na.rm=T)/(sum(res$trc_case,na.rm=T)+sum(res$trc_control,na.rm=T)))
  burden=rbind(burden,data.frame(cre=cre,type="contraction",p=p,RR=RR))  
}

get_TRGT_info_fast <- function(TRID, sampleid, hapid, motifid) {
  cmd <- sprintf(
    "bcftools view -i 'INFO/TRID=\"%s\"' %s.vcf.gz | \
     bcftools query -f '[%%SD\\t%%AP\\t%%MC\\n]' | \
     awk -v hap=%d -v motif=%d '
       {
         split($1, sd, \",\")
         split($2, ap, \",\")
         split($3, mc, \",\")
         split(mc[hap], mc_counts, \"_\")
         printf \"%%s\\t%%s\\t%%s\\n\", sd[hap], ap[hap], mc_counts[motif]
       }'",
    TRID, sampleid, hapid, motifid
  )
  
  result <- system(cmd, intern = TRUE)
  
  if (length(result) > 0) {
    values <- as.numeric(strsplit(result, "\t")[[1]])
    list(
      SD = values[1],
      AP = values[2],
      MC = values[3]
    )
  } else {
    warning("No matching record found for TRID: ", TRID)
    NULL
  }
}


exp=res[which(res$tre_case>0),c("maxtr_control","ID")]
d=d_scz[exp$ID,,drop=FALSE]
expind=data.frame()

for(i in 1:nrow(exp)){
  l=which(d[i,]>max(exp[i,1]*2,exp[i,1]+5))
  if(length(l)==0){next}
  posid=rownames(d)[i]
  pos=unlist(strsplit(posid,"_",fixed=T))
  start=as.numeric(pos[2])
  end=as.numeric(pos[3])
  chr=pos[1]
  fwrite(data.frame(chr=chr,start=start,end=end),"int.bed",sep="\t",col.names=F,scipen=20)
  system("liftOver int.bed ~/liftOver/grch38-chm13v2.chain new.bed unMapped -minMatch=0.5")
  int=data.frame()
  for (idx in l){
    TRID=exp$ID[i]
    id=colnames(d)[idx]
    id=strsplit(id,"_",fixed=T)[[1]]
    hapid=as.numeric(gsub("H","",id[2]))
    sampleid=id[1]
    motifid=as.numeric(strsplit(rownames(d_scz)[i],"_",fixed=T)[[1]][2])
    #purity and depth filtration
    info <- get_TRGT_info_fast(TRID, sampleid, hapid, motifid) 
    if (is.null(info)) { next }
    #assembly filtration
    command=paste(c("rb liftover --bed new.bed /path/",sampleid,".1.align.paf > int.paf"),collapse="")
    system(command)
    seq1="failed"
    int1=fread("int.paf",data.table=F)
    if(nrow(int1)>0){
      int1=int1[,c(1,3,4,13,14,5)]
      fwrite(int1,"int1.bed",sep="\t",col.names=F,scipen=20)
      command=paste(c("rb get-fasta --name --strand --bed int1.bed --fasta /path/",sampleid,".1.fa"),collapse="")
      seq1=system(command,intern=T)
      seq1=paste(seq1[-1],collapse="")
    }
    command=paste(c("rb liftover --bed new.bed /path/",sampleid,".2.align.paf > int.paf"),collapse="")
    system(command)
    seq2="failed"
    int1=fread("int.paf",data.table=F)
    if(nrow(int1)>0){
      int1=int1[,c(1,3,4,13,14,5)]
      fwrite(int1,"int1.bed",sep="\t",col.names=F,scipen=20)
      command=paste(c("rb get-fasta --name --strand --bed int1.bed --fasta /path/",sampleid,".2.fa"),collapse="")
      seq2=system(command,intern=T)
      seq2=paste(seq2[-1],collapse="")
    }
    int=rbind(int,data.frame(TRID=TRID,motifid=motifid,sample=sampleid,SD=info$SD,AP=info$AP,MC=info$MC,seq1=seq1,seq2=seq2))
  }
  expind=rbind(expind,int)
  print(i)
}

#candidate for IGV & assembly validation
extreme_expansion=expind[which(expind$SD>2 & expind$AP>0.7),]
