rm(list=ls())
gc()
#####network analysis##########
setwd("/Users/tq20202/Desktop/network")
info<-read.table(file="9606.protein.info.v11.5.txt", header=F, sep="\t", quote="", stringsAsFactors=F, check.names=F)
base<-read.table(file="9606.protein.links.full.v11.5.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
base<-base[base$coexpression!=0 | base$coexpression_transferred!=0,]
pool<-read.csv("pool.csv")
pool<-data.frame(gene=pool[pool$gene %in% info$V2,])
for (i in 1:nrow(pool)){
  pool[i,"id"]=info[info$V2==pool[i,"gene"],"V1"]
}
coloc<-read.csv("coloc.csv")
for (i in 1:nrow(coloc)){
  if(coloc[i,"gene"] %in% info$V2){coloc[i,"id"]=info[info$V2==coloc[i,"gene"],"V1"]}
  else{coloc[i,"id"]=NA}
}
coloc<-na.omit(coloc)
immod<-read.csv("immod.csv")
for (i in 1:nrow(immod)){
  if(immod[i,"gene"] %in% info$V2){immod[i,"id"]=info[info$V2==immod[i,"gene"],"V1"]}
  else{immod[i,"id"]=NA}
}
immod<-na.omit(immod)

#######breast cancer#############
######immod###################################################
immod<-immod[immod$type=="Immunostimulator" | immod$type=="Immunoinhibitor",]
pgene<-coloc[coloc$cancer=="prostate cancer","id"]
setz<-c()
set<-c()
for (i in 1:nrow(immod)){
  inter<-unique(c(base[base$protein1==immod[i,"id"],"protein2"],base[base$protein2==immod[i,"id"],"protein1"]))
  set<-c(set,inter)
}
set<-unique(set)
num=sum(pgene %in% set)
assoc<-0
for (k in 1:length(set)){
  t<-unique(c(base[base$protein1==set[k],"protein2"],base[base$protein2==set[k],"protein1"]))
  tt<-sum(pgene %in% t)
  assoc<-sum(assoc,tt)
}
table<-data.frame(num=num,assoc=assoc)  


n=length(pgene)
ctable<-c()
for (cicus in 1:1000){
  cgene<-sample(pool$id,size=n,replace=T)
  cset<-c()
  for (i in 1:nrow(immod)){
    inter<-unique(c(base[base$protein1==immod[i,"id"],"protein2"],base[base$protein2==immod[i,"id"],"protein1"]))
    cset<-c(cset,inter)
  }
  cset<-unique(cset)
  num=sum(cgene %in% cset)
  assoc<-0
  for (k in 1:length(cset)){
    t<-unique(c(base[base$protein1==cset[k],"protein2"],base[base$protein2==cset[k],"protein1"]))
    tt<-sum(cgene %in% t)
    assoc<-sum(assoc,tt)
  }
  temp<-data.frame(cicus=cicus,num=num,assoc=assoc)
  ctable<-rbind(ctable,temp)
}

ggplot(data=ctable, aes(x=num, y=assoc)) + geom_point()

######other MR###################################################
imcheck<-read.csv("imcheck.csv")
pgene<-coloc[coloc$cancer=="prostate cancer","id"]
setz<-c()
set<-c()
for (i in 1:nrow(imcheck)){
  inter<-unique(c(base[base$protein1==imcheck[i,"id"],"protein2"],base[base$protein2==imcheck[i,"id"],"protein1"]))
  set<-c(set,inter)
}
set<-unique(set)
num=sum(pgene %in% set)
assoc<-0
for (k in 1:length(set)){
  t<-unique(c(base[base$protein1==set[k],"protein2"],base[base$protein2==set[k],"protein1"]))
  tt<-sum(pgene %in% t)
  assoc<-sum(assoc,tt)
}
table<-data.frame(num=num,assoc=assoc)  


n=length(pgene)
ctable<-c()
for (cicus in 1:1000){
  cgene<-sample(pool$id,size=n,replace=T)
  cset<-c()
  for (i in 1:nrow(imcheck)){
    inter<-unique(c(base[base$protein1==imcheck[i,"id"],"protein2"],base[base$protein2==imcheck[i,"id"],"protein1"]))
    cset<-c(cset,inter)
  }
  cset<-unique(cset)
  num=sum(cgene %in% cset)
  assoc<-0
  for (k in 1:length(cset)){
    t<-unique(c(base[base$protein1==cset[k],"protein2"],base[base$protein2==cset[k],"protein1"]))
    tt<-sum(cgene %in% t)
    assoc<-sum(assoc,tt)
  }
  temp<-data.frame(cicus=cicus,num=num,assoc=assoc)
  ctable<-rbind(ctable,temp)
}

p<-ggplot(data=ctable, aes(x=num, y=assoc)) + geom_point()




