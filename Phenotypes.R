#!/bin/RScript
#Anna Maria Langmüller
#Cornell University
#Project: Fitness Assessment Cas9
#Purpose: Fecundity Phenotypic Assessment 

rm(list=ls())
options(stringsAsFactors = F)
my_wd<-"../../MS/Supplemental_Data/phenotype_assays/"
setwd(my_wd)

#Packages------
library(ggplot2)
library(ggpubr)
library(car)
library(lindia)
library(emmeans)
library(nortest)
library(cowplot)

#Fecundity------
d<-read.csv("raw_fecundity.csv",header=T,sep = ";")
d$eggs<-rowSums(d[,2:4],na.rm = T) #genotypes 
d$gtp<-"eGFP"
d$gtp[grep("F_",d$genotype)]<-"Cas9"
d$gtp[grep("[1|2]F/A",d$genotype)]<-"het"
d$line<-"0"
d$line[grep("1F",d$genotype)]<-"1"
d$line[grep("2F",d$genotype)]<-"2"
d$gtp<-as.factor(d$gtp)
d$line<-as.factor(d$line)
d$gtp<-relevel(d$gtp,ref="eGFP")
d$line<-relevel(d$line,ref="0")


#EDA
ggboxplot(d,x="gtp",y="eggs",color="line",add="jitter")
gghistogram(d,x="eggs",fill="gtp",alpha=0.5,add="mean")
shapiro.test(d$eggs) #normal response
leveneTest(eggs~gtp,data=d) #variance homogeneity

#MODEL
m<-lm(eggs~gtp,data = d) 
summary(m)
drop1(m,test="F") #genotype has a significant effect

plot(m) #OK
res<-resid(m)
ad.test(res) #normally distributed residuals 

cbind(orig=coef(m),confint(m)) #OK 

e<-emmeans(m,pairwise~gtp,adjust="tukey") #calculate contrasts + CI 
ci.plot<-as.data.frame(e$emmeans)
ci.plot$gtp<-factor(ci.plot$gtp,levels = c("eGFP","het","Cas9"),ordered = T)
p.values<-as.data.frame(e$contrasts)

temp<-strsplit(p.values$contrast,"-") #generate p.values data frame for plotting
p.values$group1<-sapply(temp,"[",1)
p.values$group2<-sapply(temp,"[",2)
p.values$group2<-gsub(" ","",p.values$group2)
p.values$group1<-gsub(" ","",p.values$group1)
p.values$padj<-round(p.values$p.value,3)
p.values$y.position<-c(110,120,100)
p.values$x.position<-c(1.5,2.5,2.5)
p.values

g<-ggplot(data=ci.plot,aes(x=gtp,y=emmean))+geom_bar(stat = "identity",fill="grey50",show.legend = F)+theme_classic()
g<-g+geom_pointrange(aes(x=gtp,y=emmean,ymin=lower.CL,ymax=upper.CL))
g<-g+scale_x_discrete("Female genotype",breaks=c("eGFP","het","Cas9"),labels=c("EGFP","het","Cas9_gRNAs"))
g<-g+ylab("Fecundity")+scale_y_continuous(breaks=c(0,25,50,75,100,125))
g<-g+theme(text = element_text(size=15))+coord_cartesian(ylim = c(0,130))
g<-g+stat_pvalue_manual(p.values,label = "padj",label.size = 4)
gfec<-g

cbind(coef(m),coef(m)+t(apply(X=dfbeta(m),MARGIN = 2,FUN = range))) #model is stable  

gg_cooksd(m,threshold="baseR") #non over 1 
influencePlot(m)
leveragePlot(m,"gtp") #29, 27, 1, 2, 7
outlierTest(m) #ns in Bonferroni correction

ds<-subset(d,line%in%c("1","2")) #no line effect (separate model)
ms<-lm(eggs~line,data=ds)
summary(ms)
plot(ms)
drop1(ms,test="F")


#Mate-Choice------ 
binom.test(26,38,0.5) #38 samples, 26 eGFP homozygotes offspring only 

mate<-data.frame("EGFP"=26/38,"Cas9_gRNAs"=(38-26)/38)
mate<-reshape2::melt(mate)
colnames(mate)<-c("genotype","mate")

g<-ggplot(data=mate,aes(x=genotype,y=mate))+geom_bar(stat = "identity",fill="grey50",show.legend = F)+theme_classic()
g<-g+scale_x_discrete("Male genotype",breaks=c("EGFP","Cas9_gRNAs"))
g<-g+ylab("Mate choice frequency")+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))
g<-g+theme(text = element_text(size=15))+coord_cartesian(ylim = c(0,1))
g<-g+geom_hline(yintercept = 0.5,color="black",linetype="dashed")
g<-g+annotate(geom="text",x=2,y=1,label="P = 0.033",size=4,fontface="italic")
gmate<-g

#Viability--------
via<-read.csv("raw_viability.csv")
via<-via[,1:3]
head(via)
via$r<-via$F.A/c(via$F.A+via$A.A) #calculate ratio 
head(via)
summary(via$r)
ad.test(via$r) #ratio is normally distributed 
t.test(via$r,mu=0.5,alternative = "two.sided") #ns

g<-ggplot(data=via,aes(x=r))+geom_histogram(bins = 15,color="black",fill="grey50")
g<-g+theme_classic()+coord_cartesian(xlim = c(0,1),ylim = c(0,15))
g<-g+ylab("Number of crosses")+scale_x_continuous("Fraction of heterozygous offspring",breaks = c(0,0.25,0.5,0.75,1))
g<-g+geom_vline(xintercept = 0.5,color="black",linetype="dashed")
g<-g+geom_vline(xintercept = mean(via$r),color="black",linetype="solid")
g<-g+annotate(geom="text",x=1.75/2,y=15,label="P = 0.295",size=4,fontface="italic")
g<-g+theme(text = element_text(size=15))
gvia<-g


g<-plot_grid(plotlist = list(gmate,gfec,gvia),nrow = 1,ncol = 3,labels = LETTERS[1:3])
g
ggsave("Phenotypes.png",plot = g,width = 11,height = 5,units = "in",dpi = 600)
