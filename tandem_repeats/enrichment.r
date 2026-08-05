sczfisher=data.frame()
for (pos in unique(res$pos)){
a=length(unique(large[which(large$pos==pos),"id"]))
b=length(unique(large[which(large$pos!=pos),"id"]))
c=length(unique(res[which(res$pos==pos & res$group=="none"),"id"]))
d=length(unique(res[which(res$pos!=pos & res$group=="none"),"id"]))
model=fisher.test(matrix(c(a,b,c,d),nr=2))
p=model$p.value
or=model$estimate
sczfisher=rbind(sczfisher,data.frame(pos=pos,a=a,b=b,c=c,d=d,p=p,or=or,group="large"))

a=length(unique(small[which(small$pos==pos),"id"]))
b=length(unique(small[which(small$pos!=pos),"id"]))
c=length(unique(res[which(res$pos==pos & res$group=="none"),"id"]))
d=length(unique(res[which(res$pos!=pos & res$group=="none"),"id"]))
model=fisher.test(matrix(c(a,b,c,d),nr=2))
p=model$p.value
or=model$estimate
sczfisher=rbind(sczfisher,data.frame(pos=pos,a=a,b=b,c=c,d=d,p=p,or=or,group="small"))
}

for (cre in unique(res$cre)[-1]){
a=length(unique(large[which(large$cre==cre),"id"]))
b=length(unique(large[which(large$cre!=cre),"id"]))
c=length(unique(res[which(res$cre==cre & res$group=="none"),"id"]))
d=length(unique(res[which(res$cre!=cre & res$group=="none"),"id"]))
model=fisher.test(matrix(c(a,b,c,d),nr=2))
p=model$p.value
or=model$estimate
sczfisher=rbind(sczfisher,data.frame(pos=cre,a=a,b=b,c=c,d=d,p=p,or=or,group="large"))

a=length(unique(small[which(small$cre==cre),"id"]))
b=length(unique(small[which(small$cre!=cre),"id"]))
c=length(unique(res[which(res$cre==cre & res$group=="none"),"id"]))
d=length(unique(res[which(res$cre!=cre & res$group=="none"),"id"]))
model=fisher.test(matrix(c(a,b,c,d),nr=2))
p=model$p.value
or=model$estimate
sczfisher=rbind(sczfisher,data.frame(pos=cre,a=a,b=b,c=c,d=d,p=p,or=or,group="small"))
}