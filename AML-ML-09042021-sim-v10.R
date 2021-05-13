#!/bin/RScript 
#Anna Maria Langmüller
#Cornell University
#April 2021
#Project: FAC
#Purpose: Simulate cage populations


#Packages & Functions----
rm(list=ls())
wd<-"./"
setwd(wd)
getwd()
source("AML-ML-v10.R")

library(reshape2)
library(plyr)
library(cowplot)


#Global Variables----
MA<-"data/TransitionMatrix.txt" #raw Transition Matrix path
GENOTYPES<-c("DDTT","DdTT","ddTT","DDTt","DdTt","ddTt","DDtt","Ddtt","ddtt") #Genotypes
ITER<-10000


#Empirical Data----
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

#Cas9_gRNAs----
set.seed(050689)
cm<-c(0,0,4)
om<-c(0,0,4)
ac<-rep(1,7)
e<-1
x<-1
D<-l_FACacN4
construct<-"Cas9_gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["ne"]<-175
fitness["vcs"]<-0.98
fitness["vos"]<-0.84

my_sim<-list()


for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.iter<-replicate(ITER,p0,simplify = F)
  p0.sim<-lapply(p0.iter,function(k) sim_data(p = k,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                                    fmpath = MA,ft = ac[i]))
  p0.sim.pheno<-lapply(p0.sim,function(k) t(apply(k,1,geno_to_pheno)))
  sim_res<-do.call("rbind",p0.sim.pheno)
  my_sim[[i]]<-sim_res
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,".rds",sep="")
toSave
saveRDS(my_sim,toSave)

#Cas9_no-gRNAs----
set.seed(42)
cm<-c(0,0,0)
om<-c(0,0,0)
ac<-rep(0,2)
e<-0
x<-0
D<-l_FACacN
construct<-"Cas9_no-gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["ne"]<-243

my_sim<-list()
for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.iter<-replicate(ITER,p0,simplify = F)
  p0.sim<-lapply(p0.iter,function(k) sim_data(p = k,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                                              fmpath = MA,ft = ac[i]))
  p0.sim.pheno<-lapply(p0.sim,function(k) t(apply(k,1,geno_to_pheno)))
  sim_res<-do.call("rbind",p0.sim.pheno)
  my_sim[[i]]<-sim_res
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,".rds",sep="")
toSave
saveRDS(my_sim,toSave)

#no-Cas9_no-gRNAs----
set.seed(123)
cm<-c(0,0,0)
om<-c(0,0,0)
ac<-rep(0,2)
e<-0
x<-0
D<-l_FACacR
construct<-"no-Cas9_no-gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["ne"]<-162

my_sim<-list()
for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.iter<-replicate(ITER,p0,simplify = F)
  p0.sim<-lapply(p0.iter,function(k) sim_data(p = k,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                                              fmpath = MA,ft = ac[i]))
  p0.sim.pheno<-lapply(p0.sim,function(k) t(apply(k,1,geno_to_pheno)))
  sim_res<-do.call("rbind",p0.sim.pheno)
  my_sim[[i]]<-sim_res
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,".rds",sep="")
toSave
saveRDS(my_sim,toSave)

#Cas9HF1_gRNAs----
set.seed(456)
cm<-c(0,0,0)
om<-c(0,0,0)
ac<-rep(0,2)
e<-0
x<-0
D<-l_FACacNf4
construct<-"Cas9HF1_gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["ne"]<-396

my_sim<-list()
for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.iter<-replicate(ITER,p0,simplify = F)
  p0.sim<-lapply(p0.iter,function(k) sim_data(p = k,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                                              fmpath = MA,ft = ac[i]))
  p0.sim.pheno<-lapply(p0.sim,function(k) t(apply(k,1,geno_to_pheno)))
  sim_res<-do.call("rbind",p0.sim.pheno)
  my_sim[[i]]<-sim_res
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,".rds",sep="")
toSave
saveRDS(my_sim,toSave)

#Cas9_gRNAs, full model, no drift----
set.seed(42)
cm<-c(0,0,4)
om<-c(0,0,4)
ac<-rep(1,length(l_FACacN4))
e<-1
x<-1
D<-l_FACacN4
construct<-"Cas9_gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["vcs"]<-0.98
fitness["vos"]<-0.84
fitness["ne"]<-0

my_sim<-list()
for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.sim<-sim_data(p = p0,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                                              fmpath = MA,ft = ac[i])
  p0.sim.pheno<-t(apply(p0.sim,1,geno_to_pheno))
  my_sim[[i]]<-p0.sim.pheno
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,"-no-Drift.rds",sep="")
saveRDS(my_sim,toSave)

#Cas9_gRNAs, construct model, no drift----
set.seed(42)
cm<-c(0,0,4)
om<-c(0,0,0)
ac<-rep(1,length(l_FACacN4))
e<-1
x<-1
D<-l_FACacN4
construct<-"Cas9_gRNAs"
fitness<-rep(1,9)
names(fitness)<-c("ne","vcs","vch","vos","voh","mfcs","mfch","mfos","mfoh")
fitness["vcs"]<-0.96
fitness["ne"]<-0

my_sim<-list()
for(i in c(1:length(D))){
  p0<-D[[i]][1,1:3]
  ngen<-nrow(D[[i]])-1
  p0.sim<-sim_data(p = p0,g = GENOTYPES,fne = fitness["ne"],f = fitness,fcas9 = cm,fot = om,e = e ,x = x,n = ngen,
                   fmpath = MA,ft = ac[i])
  p0.sim.pheno<-t(apply(p0.sim,1,geno_to_pheno))
  my_sim[[i]]<-p0.sim.pheno
  print(i)
}

toSave<-paste("sim",Sys.Date(),construct,paste(cm,collapse = ""),
              paste(om,collapse = ""),ac[1],x,e,sep="-")
toSave<-paste(toSave,"-no-Drift.rds",sep="")
toSave

saveRDS(my_sim,toSave)
