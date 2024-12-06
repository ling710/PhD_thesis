rm(list=ls())
gc()
setwd("/Users/tq20202/Desktop/DICE")
library(ieugwasr)
library(vcfR)
library(vroom)
library(TwoSampleMR)
library(MendelianRandomization)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

f<-list.files("data")
setwd("/Users/tq20202/Desktop/DICE/data")
data<-c()
for (i in 1:length(f)){
  temp<-read.vcfR(f[i])
  temp<-as.data.frame(temp@fix)
  temp$p=sapply(strsplit(temp$INFO,";"),"[",3)
  temp$p=as.numeric(sapply(strsplit(temp$p,"="),"[",2))
  temp$beta=sapply(strsplit(temp$INFO,";"),"[",4)
  temp$beta=as.numeric(sapply(strsplit(temp$beta,"="),"[",2))
  temp$gene=sapply(strsplit(temp$INFO,";"),"[",1)
  temp$gene=sapply(strsplit(temp$gene,"="),"[",2)
  temp$geneid=sapply(strsplit(temp$INFO,";"),"[",2)
  temp$geneid=sapply(strsplit(temp$geneid,"="),"[",2)
  temp$se=sqrt(((temp$beta)^2)/qchisq(temp$p,1,lower.tail=F))
  temp<-temp[temp$p<5e-8,]
  temp<-ld_clump_local(
    dplyr::tibble(rsid=temp$ID, chr=temp$CHROM,pos=temp$POS,beta=as.numeric(temp$beta),se=as.numeric(temp$se),effect_allele = temp$ALT,other_allele = temp$REF,pval=as.numeric(temp$p), id=temp$gene,samplesize = 91,geneid=temp$geneid), 
    clump_kb=10000,
    clump_r2=0.001,
    clump_p=1,
    bfile="/Users/tq20202/Desktop/ref/1000_genome_GRCH37/EUR_phase3_autosomes", 
    plink_bin = genetics.binaRies::get_plink_binary())
  temp$celltype=strsplit(f[i],"[.]")[[1]][1]
  data<-rbind(data,temp)
}
fre<-vroom("/Users/tq20202/Desktop/ref/fileFrequency.frq")
for (i in 1:nrow(data)){
  temp=fre[fre$ID==data[i,"SNP"],]
  if (data[i,"other_allele.exposure"]==temp$REF){data[i,"eaf"]=temp$ALT_FREQS}
  else if (data[i,"effect_allele.exposure"]==temp$ALT){data[i,"eaf"]=temp$ALT_FREQS}
  else {data[i,"eaf"]=1-temp$ALT_FREQS}
}


setwd("/Users/tq20202/Desktop/DICE")
data<-vroom("exposure.csv")
exposure_dat <- cbind(data,fstatistics=1)
for (s in 1:nrow(exposure_dat)){
  z <- exposure_dat[s,"beta.exposure"]/exposure_dat[s,"se.exposure"]
  pve <- z^2/(exposure_dat[s,"samplesize.exposure"]+z^2)
  exposure_dat[s,"fstatistics"] <- (exposure_dat[s,"samplesize.exposure"]-2)*pve/(1-pve)
}
print(min(exposure_dat$fstatistics))
print(max(exposure_dat$fstatistics))
exposure_dat<-exposure_dat[exposure_dat$fstatistics>10,]

outlist<-vroom("outlist.csv")
for (i in 1:nrow(outlist)){
  outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,outcomes = outlist[i,"IEU_GWAS_id"]$IEU_GWAS_id)
  mu=outlist[i,"N_case"]$N_case/outlist[i,"Samplesize"]$Samplesize
  if (sum(i==12:13)>0){outcome_dat$beta.outcome=outcome_dat$beta.outcome/mu;outcome_dat$se.outcome=outcome_dat$se.outcome/mu}
  harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  test<-harmdat
  test$new<- paste(test$SNP,test$exposure,sep="_")
  test<-test[,-which(colnames(test)=="id.exposure")]
  test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                      beta_col = "beta.exposure", se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",
                      other_allele_col = "other_allele.exposure",
                      eaf_col = "eaf.exposure",
                      phenotype_col = "new",samplesize_col = "samplesize.exposure")
  testout  <- outcome_dat
  inter<-intersect(test$SNP,testout$SNP)
  test <- test[test$SNP %in% inter,]
  testout <- testout[testout$SNP %in% inter,]
  testharm <- harmonise_data(test, testout, action=2)
  testharm$mr_keep<-TRUE
  testres<-mr(testharm)
  testharm<-steiger_filtering(testharm)
  unique(testharm$steiger_dir)
  
  harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
  harmdat$steiger_dir <- 0
  harmdat$steiger_pval <- 0
  for (k in 1:nrow(harmdat)){
    harmdat[k,"steiger_dir"] <- testharm[testharm$exposure==harmdat[k,"match"],"steiger_dir"][1]
    harmdat[k,"steiger_pval"] <- testharm[testharm$exposure==harmdat[k,"match"],"steiger_pval"][1]
  }
  #harmdat<-harmdat[harmdat$steiger_dir==1,]
  outtype=outlist[i,"cancer_type"]
  harmname<-paste0("/Users/tq20202/Desktop/DICE/mrharmdata/harmdata_",outtype,".txt",sep="")
  write.table(harmdat, file=harmname, row.names=F, col.names=T, sep="\t", quote=F)
  
  exp<-unique(harmdat$exposure)
  exp1=exposure_dat
  exp2=NULL
  mrres <- c()
  for (t in 1:length(exp)){
    dat <- harmdat[harmdat$exposure==exp[t],]
    dat<-dat[!is.na(dat$SNP),]
    if (nrow(dat)==1){
      result <- mr_wald_ratio(b_exp= dat$beta.exposure, b_out=dat$beta.outcome, se_exp= dat$se.exposure, se_out= dat$se.outcome)
      result <- data.frame(id=exp[t],method="Wald_ratio",nsnp=result$nsnp,b=result$b,se=result$se,CIlower=NA,CIupper=NA,
                           pval=result$pval,intercept=NA,intercept_se=NA,inter_CIlower=NA,inter_CIupper=NA,
                           intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)
    } else if (nrow(dat)==2) {
      result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                        bxse = dat$se.exposure,
                                                        by = dat$beta.outcome,
                                                        byse = dat$se.outcome))
      result <- data.frame(id=exp[t],
                           method="IVW",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError,
                           CIlower=result$CILower,
                           CIupper=result$CIUpper,
                           pval=result$Pvalue,
                           intercept=NA,
                           intercept_se=NA,
                           inter_CIlower=NA,
                           inter_CIupper=NA,
                           intercept_pval=NA,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])
    }else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){
      rho <- ld_matrix_local(variants=dat$SNP,bfile="/Users/tq20202/Desktop/ref/EUR/EUR", plink_bin = genetics.binaRies::get_plink_binary())
      result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                          bxse = dat$se.exposure,
                                                          by = dat$beta.outcome,
                                                          byse = dat$se.outcome,
                                                          cor = rho))
      result <- data.frame(id=exp[t],
                           method="MR-Egger",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError.Est,
                           CIlower=result$CILower.Est,
                           CIupper=result$CIUpper.Est,
                           pval=result$Pvalue.Est,
                           intercept=result$Intercept,
                           intercept_se=result$StdError.Int,
                           inter_CIlower=result$CILower.Int,
                           inter_CIupper=result$CIUpper.Int,
                           intercept_pval=result$Pvalue.Int,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])
      result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                         bxse = dat$se.exposure,
                                                         by = dat$beta.outcome,
                                                         byse = dat$se.outcome,
                                                         cor = rho))
      result2 <- data.frame(id=exp[t],
                            method="IVW",
                            nsnp=result2$SNPs,
                            b=result2$Estimate,
                            se=result2$StdError,
                            CIlower=result2$CILower,
                            CIupper=result2$CIUpper,
                            pval=result2$Pvalue,
                            intercept=NA,
                            intercept_se=NA,
                            inter_CIlower=NA,
                            inter_CIupper=NA,
                            intercept_pval=NA,
                            hetero_Q=result2$Heter.Stat[1],
                            hetero_pvale=result2$Heter.Stat[2])
      result <- rbind(result,result2)
    } else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){
      result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                          bxse = dat$se.exposure,
                                                          by = dat$beta.outcome,
                                                          byse = dat$se.outcome))
      result <- data.frame(id=exp[t],
                           method="MR-Egger",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError.Est,
                           CIlower=result$CILower.Est,
                           CIupper=result$CIUpper.Est,
                           pval=result$Pvalue.Est,
                           intercept=result$Intercept,
                           intercept_se=result$StdError.Int,
                           inter_CIlower=result$CILower.Int,
                           inter_CIupper=result$CIUpper.Int,
                           intercept_pval=result$Pvalue.Int,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])
      result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                         bxse = dat$se.exposure,
                                                         by = dat$beta.outcome,
                                                         byse = dat$se.outcome))
      result2 <- data.frame(id=exp[t],
                            method="IVW",
                            nsnp=result2$SNPs,
                            b=result2$Estimate,
                            se=result2$StdError,
                            CIlower=result2$CILower,
                            CIupper=result2$CIUpper,
                            pval=result2$Pvalue,
                            intercept=NA,
                            intercept_se=NA,
                            inter_CIlower=NA,
                            inter_CIupper=NA,
                            intercept_pval=NA,
                            hetero_Q=result2$Heter.Stat[1],
                            hetero_pvale=result2$Heter.Stat[2])
      result <- rbind(result,result2)}
    mrres <- rbind(mrres,result)
  }
  #####FDR Adjustment for Pvalue######
  mrres<-generate_odds_ratios(mrres)
  mrres$fdr <- p.adjust(mrres$pval, method = "fdr", n = length(mrres$pval))
  resname<-paste0("/Users/tq20202/Desktop/DICE/mr/mrres_",outtype,".txt",sep="")
  write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
}

f=list.files("mr")
setwd("./mr")
mrres<-c()
for (i in 1:length(f)){
  temp=vroom(f[i])
  temp$outcome=strsplit(strsplit(f[i],split="_")[[1]][2],split="[.]")[[1]][1]
  temp<-temp[temp$fdr<0.05,]
  mrres<-rbind(mrres,temp)
}
write.csv(mrres,"/Users/tq20202/Desktop/DICE/mrres_fdr.csv")

mrres$match<-paste(mrres$id,mrres$outcome,sep="_")
setwd("/Users/tq20202/Desktop/DICE")
f=list.files("mrharmdata")
setwd("/Users/tq20202/Desktop/DICE/mrharmdata")
colocdata<-c()
for (i in 1:length(f)){
  temp<-vroom(f[i])
  out<-strsplit(strsplit(f[i],split="_")[[1]][2],split="[.]")[[1]][1]
  temp$match<-paste(temp$exposure,out,sep="_")
  temp<-temp[temp$match %in% mrres$match,]
  colocdata<-rbind(colocdata,temp)
}
write.csv(colocdata,"colocdata.csv")

setwd("/Users/tq20202/Desktop/DICE")
colocdata<-vroom("colocdata.csv")
library(coloc)
colocres<-c()
setwd("/Users/tq20202/Desktop/DICE/data")
for (i in 364:nrow(colocdata)){
  exp=read.vcfR(colocdata[i,"file"]$file)
  exp<-as.data.frame(exp@fix)
  exp$POS<-as.numeric(exp$POS)
  snp=colocdata[i,"SNP"]$SNP
  chr=colocdata[i,"chr.exposure"]$chr.exposure
  pos=colocdata[i,"pos.exposure"]$pos.exposure
  gene=colocdata[i,"gene"]$gene
  exp<-exp[exp$CHROM==chr & exp$POS>pos-500000 & exp$POS<pos+500000,]
  exp$p=sapply(strsplit(exp$INFO,";"),"[",3)
  exp$p=as.numeric(sapply(strsplit(exp$p,"="),"[",2))
  exp$beta=sapply(strsplit(exp$INFO,";"),"[",4)
  exp$beta=as.numeric(sapply(strsplit(exp$beta,"="),"[",2))
  exp$gene=sapply(strsplit(exp$INFO,";"),"[",1)
  exp$gene=sapply(strsplit(exp$gene,"="),"[",2)
  exp<-exp[exp$gene==gene,]
  exp$se=sqrt(((exp$beta)^2)/qchisq(exp$p,1,lower.tail=F))
  exp<-exp[,-c(6:8)]
  for (k in 1:nrow(exp)){
    temp=fre[fre$ID==exp[k,"ID"],]
    if (nrow(temp)==0){next;}
    if (exp[k,"REF"]==temp$REF){exp[k,"eaf"]=temp$ALT_FREQS}
    else if (exp[k,"ALT"]==temp$ALT){exp[k,"eaf"]=temp$ALT_FREQS}
    else {exp[k,"eaf"]=1-temp$ALT_FREQS}
  }
  exp<-na.omit(exp)
  outcomeid<-colocdata[i,"id.outcome"]$id.outcome
  out<-extract_outcome_data(snps=exp$ID,outcomes = outcomeid)
  temp=intersect(exp$ID,out$SNP)
  exp<-exp[exp$ID %in% temp,]
  exp <- exp[order(exp$ID),]
  exp <- exp[!duplicated(exp$ID),]
  out<-out[out$SNP %in% temp,]
  out <- out[order(out$SNP),]
  out <- out[!duplicated(out$SNP),]
  samplesize=outlist[outlist$IEU_GWAS_id==outcomeid,"Samplesize"]$Samplesize
  ncase=outlist[outlist$IEU_GWAS_id==outcomeid,"N_case"]$N_case
  out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele.outcome,other_allele=out$other_allele.outcome,effect_allele_freq=out$eaf.outcome,beta=out$beta.outcome,se=out$se.outcome,p=out$pval.outcome,n=samplesize,case=ncase)
  exp=data.frame(SNP=exp$ID,effect_allele=exp$ALT,other_allele=exp$REF,effect_allele_freq=exp$eaf,beta=exp$beta,se=exp$se,p=exp$p,n=91)
  df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
  df <- na.omit(df)
  result <- coloc.analysis(df$SNP,df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 
  colocres<-rbind(colocres,result)
}
write.csv(colocres,"colocres.csv")