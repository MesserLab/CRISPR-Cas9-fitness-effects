#!/bin/RScript
#Anna Maria Langmueller
#for: Messer lab, Cornell University
#Project: FAC
#Purpose: empirical data, ML method 
#04/09/2020


#Functions & Packages----
rm(list=ls())
options(stringsAsFactors = F)
library(R.utils)

wd<-"./"
setwd(wd)
getwd()
source("AML-ML-v10.R") #load functions

#fD ... data
#fcm ... cas9 mode 
#fom ... off target mode 
#fm ... path raw transition matrix
#fg ... genotype names
#fac ... already cut 
#fx ... germline cut rate
#fe ... embryo cut rate
#fc_limit ... upper limit for Cas9 parameter
#fo_limit ... upper limit of off-target parameter
calc_ML<-function(fD,fcm,fom,fm,fg,fac,fmtimeout,fx,fe,fc_limit,fo_limit){
  fml<-NULL
  fml_res<-NULL
  fci<-NULL
  fci_res<-NULL
  fn<-sum(unlist(lapply(fD,function(x) nrow(x)-1))) #calculate data points 
  fp<-determine_para_range(fcm,fom,fcm_limit = fc_limit,fom_limit = fo_limit) #determine parameters
  print(fp)
  tryCatch({
    fml<-withTimeout(optim(fp[[1]],logL,method = "L-BFGS-B",control = list(fnscale=-1),lower = fp[[3]],upper = fp[[2]],fN=fD,fcm=fcm,fom=fom,ft=fac,fmpath=fm,fmyg=fg,fx=fx,fe=fe),timeout=fmtimeout)}, TimeoutException=function(ex) {
      message("Timeout in ML.Skipping.")
    })
  
  if(!is.null(fml)){
    faicc<-calculate_AICc(fml$value,length(fp[[1]]),fn) #AICc
    fml_res<-data.frame(cm=paste(fcm,collapse = ""),om=paste(fom,collapse=""),aicc=faicc,logL=fml$value)
    for(j in c(1:length(fml$par))) fml_res<-cbind(fml_res,fml$par[j])
    colnames(fml_res)[5:ncol(fml_res)]<-names(fp[[1]])
    fml_res$ac<-paste(fac,collapse="-")
    
    tryCatch({
      fci<-withTimeout(calculate_CI(fml = fml,fdist = 1/1000,falpha = 0.05,fN = fD,fcas9_modes = fcm,fot_modes = fom,g = fg,t = fac,m = fm,x = fx,e = fe),timeout=fmtimeout)}, TimeoutException=function(ex){
        message("Timeout in CI. Skipping")
      })
    if(!is.null(fci))
    {
      fci_res<-cbind(paste(fcm,collapse=""),paste(fom,collapse = ""),fci)
      colnames(fci_res)<-c("cm","om","par","low_ci","up_ci") 
    }
  }
  return(list(fml_res,fci_res))
}

#empirical data ----
#data tables (1/replicate) with counts (rr,rg,gg ) + date in fourth column
emp<-list.files(path = "data/rawData/",pattern = "*txt",full.names = T)
FACacN4<-emp[grep("FACacN4",emp)]
FACacN<-emp[grep("FACacN.txt",emp)]
FACacNf4<-emp[grep("FACacNf4",emp)]
FACacR<-emp[grep("FACacR",emp)]

l_FACacN4<-lapply(FACacN4,function(x) read.table(x,header = T))
l_FACacN<-lapply(FACacN,function(x) read.table(x,header = T))
l_FACacNf4<-lapply(FACacNf4,function(x) read.table(x,header = T))
l_FACacR<-lapply(FACacR,function(x) read.table(x,header = T))

#remove maternal effects (1FACacN4-[1|2], 2FACacN4, 5FACacN4,3FACacN4)
l_FACacN4[[1]]<-l_FACacN4[[1]][-c(1),]
l_FACacN4[[2]]<-l_FACacN4[[2]][-c(1),]
l_FACacN4[[5]]<-l_FACacN4[[5]][-c(1),]
l_FACacN4[[6]]<-l_FACacN4[[6]][-c(1),]
l_FACacN4[[7]]<-l_FACacN4[[7]][-c(1),]


#global variables----
MA<-"data/TransitionMatrix.txt" #raw Transition Matrix path
GENOTYPE<-c("DDTT","DdTT","ddTT","DDTt","DdTt","ddTt","DDtt","Ddtt","ddtt") #Genotypes
FITNESS<-"data/FitnessModes.txt"
TIMEOUT<-3000
FCLIMIT<-1
FOTLIMIT<-2

#Test option----
testing<-T
if(testing){
  print(calc_ML(fD = l_FACacN,fcm = c(0,0,0),fom = c(0,0,0),fm = MA,fac = rep(1,2),fx=1,fe=1,fg = GENOTYPE,fmtimeout=120,fc_limit = 2,fo_limit = 2))
  print(calc_ML(fD = l_FACacN,fcm = c(0,0,4),fom = c(0,0,0),fm= MA, fac= rep(1,2),fg = GENOTYPE,fmtimeout = 120,fx=0,fe=0,fc_limit = 2,fo_limit = 2))
}

#Prep Fitness modes ---- 
fitness<-read.table(FITNESS,header = F,sep = "\t")
print(fitness)

#Cas9_gRNAs, all cut ---- 
x<-1
e<-1
D<-l_FACacN4
ac<-rep(1,7)
construct<-"Cas9_gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#Cas9_gRNAs, initial off-target ----
x<-0
e<-0
D<-l_FACacN4
ac<-rep(1,7)
construct<-"Cas9_gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#Cas9_no-gRNAs, none cut ---- 
x<-0
e<-0
D<-l_FACacN
ac<-rep(0,2)
construct<-"Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#Cas9_no-gRNAs, all cut 
#---------
x<-1
e<-1
D<-l_FACacN
ac<-rep(1,2)
construct<-"Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}


#Cas9_no-gRNAs, initial off-target
#---------
x<-0
e<-0
D<-l_FACacN
ac<-rep(1,2)
construct<-"Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#no-Cas9_no-gRNAs, none cut ---- 
x<-0
e<-0
D<-l_FACacR
ac<-rep(0,2)
construct<-"no-Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#no-Cas9_no-gRNAs, all cut 
#---------
x<-1
e<-1
D<-l_FACacR
ac<-rep(1,2)
construct<-"no-Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}


#no-Cas9_no-gRNAs, initial off-target
#---------
x<-0
e<-0
D<-l_FACacR
ac<-rep(1,2)
construct<-"no-Cas9_no-gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#Cas9HF1_gRNAs, none cut ---- 
x<-0
e<-0
D<-l_FACacNf4
ac<-rep(0,2)
construct<-"Cas9HF1_gRNAs"

for(i in c(1:nrow(fitness))){
  i<-1
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

#Cas9HF1_gRNAs, all cut 
#---------
x<-1
e<-1
D<-l_FACacNf4
ac<-rep(1,2)
construct<-"Cas9-"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}


#Cas9HF1_gRNAs, initial off-target
#---------
x<-0
e<-0
D<-l_FACacNf4
ac<-rep(1,2)
construct<-"Cas9HF1_gRNAs"

for(i in c(1:nrow(fitness))){
  cas9_mode<-as.numeric(unlist(fitness[i,c(1:3)]))
  ot_mode<-as.numeric(unlist(fitness[i,c(4:6)]))
  res<-calc_ML(fD = D, fcm = cas9_mode,fom = ot_mode,fm=MA,fg = GENOTYPE,fac = ac,fmtimeout = TIMEOUT,
               fx = x,fe = e,fc_limit = FCLIMIT,fo_limit = FOTLIMIT)
  toStore<-paste(Sys.Date(),construct,paste(cas9_mode,collapse = ""),paste(ot_mode,collapse = ""),collapse = "",ac[1],x,e,sep="-")
  toStore<-paste(toStore,".rds",sep="")
  saveRDS(res,toStore)
}

