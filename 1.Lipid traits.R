ao2<-ao%>%
  subset(grepl('24097068',pmid,ignore.case = TRUE))
ao3<-ao%>%
  subset(grepl('34226706',pmid,ignore.case = TRUE))
ao4<-ao%>%
  subset(grepl('34059833',pmid,ignore.case = TRUE))
ao5<-ao%>%
  subset(grepl('27989323',pmid,ignore.case = TRUE))
ao6<-ao%>%
  subset(grepl('26343387',pmid,ignore.case = TRUE))

out<-fread("~/R/GWASR9FINNGEN/finngen_R9_M13_OSTEONECROSIS.gz",header = T)#读取本地结局（这里将文件名换成你自己本地结局的文件名）
out$trait <- 'Osteonecrosis'
head(out)
out$n<-359399
outcomeid <- out
rm(out)
head(outcomeid)
outcomeid<-dplyr::rename(outcomeid,SNP=rsids)

outcomedata<-function(Lipid_trait,outcomeid){
  outcome_dat<-merge(Lipid_trait,outcomeid,by.x='SNP',by.y='SNP')
  write.csv(outcome_dat,file = "d.csv")
  out <-read_outcome_data(
    snps = Lipid_trait$SNP, 
    filename = "d.csv",
    sep = ",",
    samplesize_col = 'n',
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval",gene_col = 'nearest_genes',
    phenotype_col = 'trait',
    chr_col = '#chrom',
    pos_col = 'pos')
  return(out)
}


# 1.Lipid traits and Rotator cuff syndrome risk####

#脂质指标
LDLC<-extract_instruments(outcome="ieu-a-300")
TG<-extract_instruments(outcome="ieu-a-302")
TC<-extract_instruments(outcome="ieu-a-301")
HDL<-extract_instruments(outcomes='ieu-a-299')
ApoA1<-extract_instruments(outcomes='ebi-a-GCST90025955')
ApoB<-extract_instruments(outcomes='ebi-a-GCST90025952')
Lpa<-extract_instruments(outcomes='ebi-a-GCST90025993')

#全身疾病
T2DM<-extract_instruments(outcome='ebi-a-GCST006867')
glucocorticoids<-extract_instruments('ebi-a-GCST90019000')
SLE<-extract_instruments('ebi-a-GCST003156')


#炎症指标
CRP<-extract_instruments(outcome= 'ebi-a-GCST90029070')

#血糖指标
Twohour_glucose<-extract_instruments(outcome= 'ebi-a-GCST90002227')
Fasting_insulin<-extract_instruments(outcome= 'ebi-a-GCST90002238')
Fasting_glucose<-extract_instruments(outcome= 'ebi-a-GCST90002232')
Glycated_hemoglobin<-extract_instruments(outcome= 'ebi-a-GCST90002244')
UA<-extract_instruments(outcome='ebi-a-GCST90018977')

#血压指标
DBP<-extract_instruments(outcome='ebi-a-GCST90000059')
Pulse_pressure<-extract_instruments(outcome='ebi-a-GCST90000061')
SBP<-extract_instruments(outcome='ebi-a-GCST90000062')


outstrings<-paste0('ebi-a-GCST0044',c(20:60))
label<-ao%>%
  subset(grepl('27989323',pmid,ignore.case = TRUE))
expclump<- extract_instruments(outstrings,p1 = 5e-06,
                               clump = TRUE,
                               p2 = 5e-06,
                               r2 = 0.001,
                               kb = 10000)
expclump$exposure[expclump$id.exposure %in% label$id] <- label$trait[match(expclump$id.exposure, label$id)]
table(expclump$exposure)
cytokines<-expclump
table(cytokines$exposure)

combine1 <- do.call(rbind, lapply(list(LDLC, TG, TC, HDL, ApoA1, ApoB, Lpa), as.data.frame))
combine2 <- do.call(rbind, lapply(list(T2DM,glucocorticoids,SLE,CRP,
                                       Twohour_glucose,Fasting_insulin,Fasting_glucose,
                                       Glycated_hemoglobin,UA,DBP,Pulse_pressure,SBP), as.data.frame))
combine3<-do.call(rbind, lapply(list(combine2,cytokines), as.data.frame))

combine<-rbind(combine1,combine2)
table(combine$exposure)



Lipid_trait<-combine1
Lipid_trait <- Lipid_trait[complete.cases(Lipid_trait),]

library(FastTraitR)
snppheno<-look_trait(Lipid_trait$SNP, pval=5e-8)
table(snppheno$trait)
confounder_out <- function(exp, snppheno1, cf) {
  patterns <- paste(cf, collapse = "|")
  pp <- grep(patterns, snppheno$trait, ignore.case = TRUE)
  abc <- exp[!(exp$SNP %in% snppheno$rsid[pp]), ]
  return(abc)
}
Lipid_trait<-confounder_out(Lipid_trait,snppheno,c('smoking','body mass index','alcohol','diabetes'))
snpout<-setdiff(combine1$SNP,Lipid_trait$SNP)
table(Lipid_trait$exposure)

table(Lipid_trait$exposure)
Lipid_trait<-get_f(Lipid_trait)
