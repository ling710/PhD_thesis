rm(list=ls())
gc()
library(TwoSampleMR)
library(MendelianRandomization)
library(vroom)
setwd("/Users/tq20202/Desktop/oneK1K")
fname<-list.files("data")
setwd("/Users/tq20202/Desktop/oneK1K/data")
data<-c()
for (i in 1:length(fname)){
  exposure_dat<-vroom(fname[i])
  exposure_dat$A1_FREQ<-1-exposure_dat$A2_FREQ_ONEK1K
  exposure_dat$beta<-exposure_dat$SPEARMANS_RHO/sqrt(2*exposure_dat$A1_FREQ*(1-exposure_dat$A1_FREQ))
  exposure_dat$z<-qnorm(exposure_dat$P_VALUE,lower.tail = F)
  exposure_dat$se <- abs(exposure_dat$beta/exposure_dat$z)
  exposure_dat$se=exposure_dat$se+0.00001
  exposure_dat<-exposure_dat[exposure_dat$P_VALUE<=5e-8,]
  exposure_dat$pve<-exposure_dat$z^2/(1267758+exposure_dat$z^2)
  exposure_dat$fstatistics<-(1267758-2)*exposure_dat$pve/(1-exposure_dat$pve)
  exposure_dat<-exposure_dat[exposure_dat$fstatistics>10,]
  tmp <- data.frame(SNP = exposure_dat$RSID,beta = exposure_dat$beta,se = exposure_dat$se,
                    effect_allele = exposure_dat$A1,eaf = exposure_dat$A1_FREQ,
                    other_allele = exposure_dat$A2,pval = exposure_dat$P_VALUE,
                    chr = exposure_dat$CHR,pos = exposure_dat$POS,samplesize = 1267758,
                    phenotype_id = paste(exposure_dat$GENE,exposure_dat$CELL_TYPE,sep="_"))
  exposure_dat <- format_data(tmp, type="exposure",effect_allele_col = "effect_allele",other_allele_col = "other_allele",phenotype_col = "phenotype_id",samplesize_col = "samplesize",chr_col = "chr",pos_col = "pos")
  exposure_dat <-clump_data(exposure_dat,clump_kb=10000,clump_r2=0.001,clump_p1=1,clump_p2=1,pop="EUR")
  data<-rbind(data,exposure_dat)
} 
write.csv(data,"/Users/tq20202/Desktop/oneK1K/exposure.csv")

exposure_dat<-read.csv("/Users/tq20202/Desktop/oneK1K/exposure.csv")
setwd("/Users/tq20202/Desktop/oneK1K")
outcomelist<-read.csv("outcome_selection.csv")
outcome_id <-outcomelist$IEU.GWAS.id
outcome_name<-outcomelist$cancer.type
outcome_samplesize<-outcomelist$Sample.size
ncase<-outcomelist$N.case
ncontrol<-outcomelist$N.control
for (j in 1:length(outcome_id)){
  outcome <- extract_outcome_data(snps=exposure_dat$SNP, outcomes= outcome_id[j])
  outcome<-outcome[!duplicated(outcome$SNP),]
  inter<-intersect(exposure_dat$SNP,outcome$SNP)
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
  outcome_dat <- outcome[outcome$SNP %in% inter,]
  harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  #harmdat$mr_keep<-TRUE  #consider ambiguous
  
  test<-harmdat
  test$new<- paste(test$SNP,test$exposure,sep="_")
  test<-test[,-which(colnames(test)=="id.exposure")]
  test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                      beta_col = "beta.exposure", se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",
                      eaf_col = "eaf.exposure",
                      other_allele_col = "other_allele.exposure",
                      phenotype_col = "new",samplesize_col = "samplesize.exposure")
  testout  <- outcome_dat
  inter<-intersect(test$SNP,testout$SNP)
  test <- test[test$SNP %in% inter,]
  testout <- testout[testout$SNP %in% inter,]
  testharm <- harmonise_data(test, testout, action=2)
  testharm$mr_keep<-TRUE
  testres<-mr(testharm)
  ###steiger filtering
  testres <- generate_odds_ratios(testres)
  testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.out=1,p.exp=1,p.out=1,samplesize.exposure=1)
  for (r in 1: nrow(testres)){
    testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
    testres[r,(ncol(testres)-7):ncol(testres)] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                  "se.outcome","eaf.outcome",
                                                                                  "pval.exposure","pval.outcome","samplesize.exposure")]
  }
  #testres <- na.omit(testres)
  testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.outcome=outcome_samplesize[j])
  testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=testres$samplesize.exposure))^2
  testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase[j],ncontrol=ncontrol[j],prevalence=ncase[j]/outcome_samplesize[j]))^2
  st <- psych::r.test(
    n = testres$samplesize.exposure, 
    n2 = testres$samplesize.outcome, 
    r12 = sqrt(testres$rsq.exposure), 
    r34 = sqrt(testres$rsq.outcome)
  )
  testres$steiger_dir <- testres$rsq.exposure > testres$rsq.outcome
  testres$steiger_pval <- pnorm(-abs(st$z)) * 2
  
  harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
  harmdat$steiger_dir <- 0
  harmdat$steiger_pval <- 0
  for (k in 1:nrow(harmdat)){
    harmdat[k,"steiger_dir"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_dir"][1]
    harmdat[k,"steiger_pval"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_pval"][1]
  }
  harmdat<-harmdat[harmdat$steiger_dir==1,]
  harmdat<-harmdat[harmdat$mr_keep==T,]
  length(unique(harmdat$exposure))
  harmname<-paste("harmdata",outcome_name[j],outcome_id[j],sep="_")
  harmname<-paste(harmname,".txt",sep="")
  write.table(harmdat, file=harmname, row.names=F, col.names=T, sep="\t", quote=F)
  
  #Generalized IVW  and MR-Egger for the "MendelianRandomization" package
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
                           intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)}
    else if (nrow(dat)==2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
    result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                      bxse = dat$se.exposure,
                                                      by = dat$beta.outcome,
                                                      byse = dat$se.outcome,
                                                      cor = rho))
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
                         hetero_pvale=result$Heter.Stat[2])}
    else if (nrow(dat)==2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
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
                         hetero_pvale=result$Heter.Stat[2])}
    else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
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
    result <- rbind(result,result2)}
    else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
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
  resname<-paste("mrres",outcome_name[j],outcome_id[j],sep="_")
  resname<-paste(resname,".txt",sep="")
  write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
}

#####FDR result########
setwd("/Users/tq20202/Desktop/oneK1K")
f<-list.files("mrres")
setwd("/Users/tq20202/Desktop/oneK1K/mrres")
data<-c()
for (i in 1:length(f)){
  temp<-read.table(file=f[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$cancer<-strsplit(f[i],split="_")[[1]][2]
  temp$cancerid<-strsplit(strsplit(f[i],split="cancer_")[[1]][2],split="[.]")[[1]][1]
  data<-rbind(data,temp)
}
res<-data[data$fdr<0.05,]
write.csv(res,"/Users/tq20202/Desktop/oneK1K/mrresfdr.csv")
res$match<-paste(res$id,res$cancer,res$cancerid,sep="_")
setwd("/Users/tq20202/Desktop/oneK1K")
f<-list.files("harmdata")
setwd("/Users/tq20202/Desktop/oneK1K/harmdata")
hdata<-c()
for (i in 1:length(f)){
  temp<-read.table(file=f[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$cancer<-strsplit(f[i],split="_")[[1]][2]
  temp$cancerid<-strsplit(strsplit(f[i],split="cancer_")[[1]][2],split="[.]")[[1]][1]
  hdata<-rbind(hdata,temp)
}
hdata$match<-paste(hdata$exposure,hdata$cancer,hdata$cancerid,sep="_")
adata<-hdata[hdata$match %in% res$match,]
adata$gene<-sapply(strsplit(as.character(adata$exposure),'_'), "[", 1)
adata$cell<-sapply(strsplit(as.character(adata$exposure),'_'), "[", 2)
colocdata<-data.frame(SNP=adata$SNP,chr=adata$chr,pos=adata$pos,beta.exp=adata$beta.exposure,
                      se.exp=adata$se.exposure,p.exp=adata$pval.exposure,alt.exp=adata$effect_allele.exposure,
                      ref.exp=adata$other_allele.exposure,gene=adata$gene,cell=adata$cell,
                      cancer=adata$cancer,cancerid=adata$cancerid)
write.csv(colocdata,"/Users/tq20202/Desktop/oneK1K/colocdata.csv")


##################coloc###################################
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(vroom)
setwd("/Users/tq20202/Desktop/oneK1K")

##coloc function##
coloc.analysis <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s){
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  #type cc (case control) for binary study
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="cc",s=s,N=N2)
  #setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12)
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  return(df)
}

res<-read.csv("colocdata.csv")
outcomelist<-read.csv("outcome_selection.csv")
outcome_id <-outcomelist$IEU.GWAS.id
outcome_name<-outcomelist$cancer.type
outcome_samplesize<-outcomelist$Sample.size
ncase<-outcomelist$N.case
ncontrol<-outcomelist$N.control
bres<-c()
for (i in 1: nrow(res)){
  snp<-res[i,"SNP"]
  chr<-res[i,"chr"]
  pos<-res[i,"pos"]
  cell<-res[i,"cell"]
  gene=res[i,"gene"]
  cellfile<-paste("/Users/tq20202/Desktop/oneK1K/data/",res[i,"cellfile"],sep="")
  cancer<-res[i,"cancer"]
  outid<-res[i,"cancerid"]
  f<-vroom(file=cellfile)
  f<-f[f$CHR==chr & f$POS>pos-500000 & f$POS<pos+500000,]
  f$A1_FREQ<-1-f$A2_FREQ_ONEK1K
  f$beta<-f$SPEARMANS_RHO/sqrt(2*f$A1_FREQ*(1-f$A1_FREQ))
  f$z<-qnorm(f$P_VALUE,lower.tail = F)
  f$se <- abs(f$beta/f$z)
  f$se<-f$se+0.00001
  f$n=1267758
  exp=data.frame(SNP=f$RSID,effect_allele=f$A1,effect_allele_freq=f$A1_FREQ,other_allele=f$A2,beta=f$beta,se=f$se,p=f$P_VALUE,n=f$n)
  out <- extract_outcome_data(snps=exp$SNP, outcomes= outid)
  temp=intersect(exp$SNP,out$SNP)
  exp<-exp[exp$SNP %in% temp,]
  exp <- exp[order(exp$SNP),]
  exp <- exp[!duplicated(exp$SNP),]
  out<-out[out$SNP %in% temp,]
  out <- out[order(out$SNP),]
  out <- out[!duplicated(out$SNP),]
  out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele.outcome,other_allele=out$other_allele.outcome,effect_allele_freq=out$eaf.outcome,beta=out$beta.outcome,se=out$se.outcome,p=out$pval.outcome,n=outcome_samplesize[which(outcome_id==outid)],case=ncase[which(outcome_id==outid)])
  expfile<-paste("/Users/tq20202/Desktop/oneK1K/pwcoco/",i,"_exp.txt",sep="")
  outfile<-paste("/Users/tq20202/Desktop/oneK1K/pwcoco/",i,"_out.txt",sep="")
  write.table(exp, file=expfile, row.names=F, col.names=T, sep="\t", quote=F)
  write.table(out, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
  df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
  df <- na.omit(df)
  #coloc analysis
  result <- coloc.analysis(df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 
  result <- data.frame(result,num=i,gene=gene,cell=cell,cancer=cancer,cancerid=outid)
  bres <- rbind(bres,result)
}
bres$match<-paste(bres$gene,bres$cell,bres$cancer,bres$cancerid,sep="_")
mr<-read.csv("mrresfdr.csv")
mr$match<-paste(mr$id,mr$cancer,mr$cancerid,sep="_")
for (i in 1:nrow(mr)){
  if (nrow(bres[bres$match==mr[i,"match"],])==1){mr[i,"PPH4"]<-bres[bres$match==mr[i,"match"],"PP.H4.abf"]}
  else {mr[i,"PPH4"]<-max(bres[bres$match==mr[i,"match"],"PP.H4.abf"])}
}
write.table(bres, file="/Users/tq20202/Desktop/oneK1K/bres.txt", row.names=F, col.names=T, sep="\t", quote=F)


#####LDCHECK#######################################
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)
setwd("/Users/tq20202/Desktop/oneK1K")
fdrres<-read.csv("colocdata.csv")
ldres<-c()
for (i in 1:nrow(fdrres)){
  rsid=fdrres[i,"SNP"]
  chr=fdrres[i,"chr"]
  pos=fdrres[i,"pos"]
  gene=fdrres[i,"gene"]
  cell=fdrres[i,"cell"]
  cellfile<-paste("/Users/tq20202/Desktop/oneK1K/data/",res[i,"cellfile"],sep="")
  cancer<-res[i,"cancer"]
  outid<-res[i,"cancerid"]
  f<-vroom(file=cellfile)
  f<-f[f$CHR==chr & f$POS>pos-500000 & f$POS<pos+500000,]
  assoc<-extract_outcome_data(snps=f$RSID, outcomes= outid)
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$pval.outcome),]
    data <- assoc[assoc$pval.outcome<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data[,1:9])
      snp <- append(data$SNP, rsid)
      a <- NULL
      attempts <- 0
      while(attempts<=10){    
        a <- ld_matrix(variants=snp,pop = "EUR")
        if(is.null(a)){attempts<-attempts+1}else{break}
      }
      if(is.null(nrow(a))==TRUE){c<-cbind(snp=rsid, ld_snp=rsid, ld_r2=1)} 
      else {col.index <- which(grepl(rsid,colnames(a)))
      if (length(col.index)>0){
        b <- (a[,col.index])^2
        b <- b[order(b)]
        b <- b[(length(b)-1)] 
        c <- cbind(snp=rsid, ld_snp=names(b), ld_r2=as.numeric(b))
      } else {c <- cbind(snp=rsid, ld_snp="NA", ld_r2="NA")}
      }
    } else {
      c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  } else {c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  
  ldres <- rbind(ldres,c)
}
ldres$match<-bres$match
for (i in 1:nrow(mr)){
  if (nrow(ldres[ldres$match==mr[i,"match"],])==1){mr[i,"ldcheck"]<-ldres[ldres$match==mr[i,"match"],"ld_r2"]}
  else {mr[i,"ldcheck"]<-max(ldres[ldres$match==mr[i,"match"],"ld_r2"])}
}
write.csv(ldres,"/Users/tq20202/Desktop/oneK1K/ldres.csv")
write.csv(mr,"/Users/tq20202/Desktop/oneK1K/mrresfdr.csv")

#############PWCOCO####################################
t<-c()
s1<-"./pwcoco --bfile /user/home/tq20202/pwcoco/pwcoco-master/build/reference_panal/EUR/EUR  --sum_stats1 oneK1K/"
s2<-"_exp.txt  --sum_stats2  oneK1K/"
s3<-"_out.txt  --out oneK1Kres/"
s4<-" --out_cond"
for(i in 1:44){
  temp<-paste(s1,i,s2,i,s3,i,s4,sep="")
  t<-rbind(t,temp)
}
write.csv(t,"/Users/tq20202/Desktop/oneK1K/pwcoco.csv")

setwd("/Users/tq20202/Desktop/oneK1K")
f=list.files("pwcocores")
setwd("/Users/tq20202/Desktop/oneK1K/pwcocores")
pwres<-c()
for (i in 1:length(f)){
  temp<-read.delim(f[i], header=T)
  temp<-temp[temp$H4==max(temp$H4),][1,]
  pwres<-rbind(pwres,temp)
}
write.csv(pwres,"/Users/tq20202/Desktop/oneK1K/pwres.csv")












                           