#Anna Maria Langmüller
#September 2020
#Messer lab, Cornell University
#ML Framework Functions


require(stringr)

#Functions------


#rho (see Genetics 2019)
#obs ... observed phenotype frequencies (rr,rg,gg)
#pred ... predicted phenotype frequencies (rr,rg,gg)
#n ... effective population size 
#return(r)... l for likelihood 
rho <- function(obs,pred,n) {
  Obs<-obs*n
  log.a<-2*log(n)-log(1+(pred[1]^n+pred[2]^n+pred[3]^n)/6-((1-obs[1])^n+(1-obs[2])^n+(obs[1]+obs[2])^n)/2)
  r <-  log.a +
    lgamma(n+1)-lgamma(Obs[1]+1)-lgamma(Obs[2]+1)-lgamma(Obs[3]+1) 
  
  if(pred[1]>0){
    r <- r+Obs[1]*log(pred[1])
  }
  if(pred[2]>0){
    r <- r+Obs[2]*log(pred[2])
  }
  if(pred[3]>0){
    r <- r+Obs[3]*log(pred[3])
  }
  return(r)
}

#init_m
#x ... raw transition matrix (counts)
#y ... names of count columns
#return(z) ... expected offspring frequencies
init_m<-function(x,y){
  z<-read.table(x,header = T)
  z<-cbind(z[,1:2],z[,c(3:11)]/rowSums(z[,c(3:11)]))
  colnames(z)<-c("m","f",y)
  return(z)
}

#init_genotype
#x ... names of genotypes
#return(y) ... empty named genotype vector
init_genotype<-function(x){
  y<-rep(0,9)
  names(y)<-x
  return(y)
}

#normalize_counts --> maybe get rid of that; normalize counts on site!
#x ... data frame with phenotype counts 
#return(y) ... matrix with phenotype freqs
normalize_counts<-function(x){
  y<-x[,1:3]
  y<-y/rowSums(y)
  colnames(y)<-c("rr","rg","gg")
  return(y)
}

#obtain_phenoID
#x ... named vector of genotypes (m+f, length=18)
#return(y) ... list of genotype->phenotype (3x6)
obtain_phenoID<-function(x){
  r<-grep("d",names(x))
  g<-grep("D",names(x))
  h<-intersect(r,g)
  r<-setdiff(r,h)
  g<-setdiff(g,h)
  y<-list(r,h,g)
  names(y)<-c("rr","rg","gg")
  return(y)
}

#germ_cut_d 
#v ... original offspring frequencies (length = 9)
#x ... cut rate T -> t in germline; one/two drive alleles act the same
#return(w) ... adapted offspring frequencies
germ_cut_d<-function(v,x){
  w<-v
  w["DdTT"]<-v["DdTT"]*(1-x)^2
  w["ddTT"]<-v["ddTT"]*(1-x)^2
  w["DdTt"]<-v["DdTt"]*(1-x)+v["DdTT"]*2*x*(1-x)
  w["ddTt"]<-v["ddTt"]*(1-x)+v["ddTT"]*2*x*(1-x)
  w["Ddtt"]<-v["DdTT"]*x^2+v["DdTt"]*x+v["Ddtt"]
  w["ddtt"]<-v["ddtt"]+v["ddTT"]*x^2+v["ddTt"]*x
  return(w)
}

#embryo_cut_d
#v ... original offspring frequencies (length = 9)
#e .. cut rate T -> t if mom had at least one drive allele 
#return(w) ... adapted offspring frequencies 
embryo_cut_d<-function(v,e){
  w<-v
  w["DDTT"]<-v["DDTT"]*(1-e)^2
  w["DdTT"]<-v["DdTT"]*(1-e)^2 
  w["ddTT"]<-v["ddTT"]*(1-e)^2
  w["DDTt"]<-v["DDTT"]*2*e*(1-e)+v["DDTt"]*(1-e)
  w["DdTt"]<-v["DdTT"]*2*e*(1-e)+v["DdTt"]*(1-e)
  w["ddTt"]<-v["ddTT"]*2*e*(1-e)+v["ddTt"]*(1-e)
  w["DDtt"]<-v["DDTT"]*e^2+v["DDTt"]*e+v["DDtt"]
  w["Ddtt"]<-v["DdTT"]*e^2+v["DdTt"]*e+v["Ddtt"]
  w["ddtt"]<-v["ddTT"]*e^2+v["ddTt"]*e+v["ddtt"]
  return(w)
}

#mod_transition_m
#m ... transition matrix 
#x ... germline cut rate 
#e ... embryo cut rate 
#g ... names of genotypes
#return(y) ... transition matrix with adapted offspring frequencies after germline + embryo cutting
mod_transition_m<-function(m,x,e,g){
  y<-m
  #Step 1: Germline cuts 
  for(i in c(1:nrow(y))){
    fy_iter<-y[i,] #store raw  transition probabilities + cross
    fmale_genotype<-init_genotype(x = g) #empty male genotype vector
    ffemale_genotype<-init_genotype(x = g) #empty female genotype vector
    fmale_genotype[fy_iter[["m"]]]<-1 #set male genotype
    ffemale_genotype[fy_iter[["f"]]]<-1 #set female genotype
    fmale_genotype<-germ_cut_d(fmale_genotype,x)  #apply germline cuts in males
    ffemale_genotype<-germ_cut_d(ffemale_genotype,x) #apply germline cuts in females
    fnew_mates<-as.vector(outer(ffemale_genotype,fmale_genotype))  #all possible mating probs 
    fnew_progeny<-fnew_mates*m[,(3:11)] #determine new weights
    y[i,(3:11)]<-colSums(fnew_progeny)/sum(fnew_progeny) #replace marginal probs for offspring
  }
  
  #Stemp 2: Embryo cuts
  fd_carrier<-grep("d",colnames(y))-2
  for(i in c(1:nrow(y))){
    if(y$f[i]%in%fd_carrier){
      y[i,3:11]<-embryo_cut_d(y[i,3:11],e)
    }
  }
  return(y)
}

#geno_to_pheno
#x ... named vector of genotypes (m+f,length=18)
#return= vector of phenotype frequencies | x (length=3)
geno_to_pheno<-function(x){
  i<-obtain_phenoID(x)
  r<-sum(x[i[["rr"]]])
  g<-sum(x[i[["gg"]]])
  h<-sum(x[i[["rg"]]])
  y<-c(r,h,g)
  y<-y/sum(y)
  names(y)<-c("rr","rg","gg")
  return(y)
}

#pheno_to_geno
#graw ... named vector of genotypes (m+f, length=18)
#pobs ... observed phenotypes (length=3)
#return(fgcorr) ... vector of genotype frequencies, corrected by observed phenotypes (m+f,length=18)
pheno_to_geno<-function(g,p){
  fpobs<-p/sum(p)
  fpraw<-geno_to_pheno(g) #marginal sums of genotypes giving the same phenotype 
  i<-obtain_phenoID(g)
  fgcorr<-numeric(length = 18)
  names(fgcorr)<-names(g)
  #standardize 
  for(j in i[["rr"]]) fgcorr[j]<-g[j]/fpraw["rr"]*fpobs["rr"]
  for(j in i[["rg"]]) fgcorr[j]<-g[j]/fpraw["rg"]*fpobs["rg"]
  for(j in i[["gg"]]) fgcorr[j]<-g[j]/fpraw["gg"]*fpobs["gg"]
  return(fgcorr)
}

#start_genotype .... function to start a genotype 
#x ... counts 
#fg ... genotype names
#t ... already cut rate
#return(y)... genotype freques males + females 
start_genotype<-function(x,fg,t){
  x<-unlist(x)
  names(x)<-c("rr","rg","gg")
  x<-x/sum(x)
  y<-numeric(length = 9)
  names(y)<-fg
  id<-obtain_phenoID(y)
  #none wt|wt are already cut
  y[id[["gg"]][1]]<-x["gg"]/2
  #cutting happened according to cut-rate
  y[id[["rg"]][1]]<-x["rg"]*(1-t)^2/2
  y[id[["rg"]][2]]<-x["rg"]*2*t*(1-t)/2
  y[id[["rg"]][3]]<-x["rg"]*t^2/2
  y[id[["rr"]][1]]<-x["rr"]*(1-t)^2/2
  y[id[["rr"]][2]]<-x["rr"]*2*t*(1-t)/2
  y[id[["rr"]][3]]<-x["rr"]*t^2/2
  y<-rep(y,2)
  return(y)
}

#determine_fitness_coefficients
#mode ... determines the fitness costs applied 
#fg ... names of genotypes
#a ... either d/t to determine allele
#p... parameter
#return(x) ... fitness vector for genotypes (n=9); rel weights
determine_fitness_coefficients<-function(m,fg,a,p){
  x<-rep(1,9)
  names(x)<-fg
  fa_count<-str_count(names(x),a)
  if(length(p)!=2) return("Not right parameter dimension")	  
  if(m==0){ #neutral -> nothing
  }
  if(m==1){ #dominant
    x[which(fa_count>=1)]<-p[1]
  } 
  if(m==2){ #co-dominant
    x[which(fa_count==2)]<-p[1]
    x[which(fa_count==1)]<-sqrt(p[1])
  }
  if(m==3){ #Aa != aa 
    x[which(fa_count==2)]<-1-p[1]
    x[which(fa_count==1)]<-1-(p[1]*p[2])
  }
  if(m==4){ #multiplicative
    x[which(fa_count==2)]<-p[1]^2
    x[which(fa_count==1)]<-p[1]
  }
  return(x)
}

#propagate_genotypes
#p ... parental genotypes (18 (m+f))
#m ... transition matrix (germ + embryo cuts + fitness)
#g ... names of genotypes
#fcas9 ... cas9 fitness cost
#fot ... off target fitness cost
#f... fitness parameter (named!)
#return=expected genotype frequencies next generation (18 (m+f))
propagate_genotypes<-function(p,m,g,fcas9,fot,f){
  fmale<-p[1:9]/sum(p[1:9]) #standardize male freqs
  ffemale<-p[10:18]/sum(p[10:18]) #standardize female freqs
  
  #mate choice selection 
  fmate_cas9<-determine_fitness_coefficients(m = fcas9[1],fg = g,a = "d",p = f[c("mfcs","mfch")])
  fmate_ot<-determine_fitness_coefficients(m = fot[1],fg = g,a = "t",p = f[c("mfos","mfoh")])
  fmate<-fmate_cas9*fmate_ot
  fmale<-fmale*fmate/sum(fmale*fmate)
  fpairs<-as.vector(outer(ffemale,fmale))
  fmmate<-fpairs*m[,3:11]
  
  #fecundity selection
  ffec_cas9<-determine_fitness_coefficients(m = fcas9[2],fg = g,a = "d",p = f[c("mfcs","mfch")])
  ffec_cas9<-rep(ffec_cas9,length(fmale))
  ffec_ot<-determine_fitness_coefficients (m = fot[2],fg = g,a = "t",p = f[c("mfos","mfoh")])
  ffec_ot<-rep(ffec_ot,length(fmale))
  ffec<-ffec_cas9*ffec_ot
  fmfec<-ffec*fmmate
  
  #viability selection
  fvia_cas9<-determine_fitness_coefficients(m = fcas9[3],fg = g,a = "d",p = f[c("vcs","vch")])
  fvia_ot<-determine_fitness_coefficients(m=fot[3],fg = g,a = "t",p = f[c("vos","voh")])
  fvia<-fvia_cas9*fvia_ot
  fmvia<-t(fvia*t(fmfec))
  
  #extend to males + females again 
  g<-colSums(fmvia)/sum(fmvia)
  g<-rep(g,2)
  g<-g/sum(g)
  return(g)

}

#sim_data
#p ... parental phenotype counts 
#g ... names of genotypes
#fne ... effective population size
#f .... fitness params
#fcas9 ... cas9 fitness cost
#fot ... off-target fitness cost
#e ... embryo cut rate
#x ... germline cut rate
#n ... number of generations to simulate
#fmpath ... path to transition matrix 
#ft ... already cut parameter
sim_data<-function(p,g,fne,f,fcas9,fot,e,x,n,fmpath,ft){
  r<-p/sum(p)
  r<-start_genotype(r,g,ft)
  #print(r)
  y<-r
  m<-init_m(x = fmpath,y = g) #raw transition matrix
  fmmod<-mod_transition_m(m = m,x = x,e = e,g = g) #germline/embryo cut
  for(i in c(1:n)){ #for each generation 
    fp1<-propagate_genotypes(p=r,m=fmmod,g=g,fcas9=fcas9,fot=fot,f=f)
    fpm<-fp1[1:9]/sum(fp1[1:9])
    fpf<-fp1[10:18]/sum(fp1[10:18])
    if(fne>0){
      gm<-t(rmultinom(1,fne/2,fpm))
      gm<-gm/sum(gm)
      gf<-t(rmultinom(1,fne/2,fpf))
      gf<-gf/sum(gf)
    } else {
      gm<-fpm
      gf<-fpf
    }
    r<-c(gm,gf)
    r<-r/sum(r)
    y<-rbind(y,r)
  }
  rownames(y)<-c(0:n)
  colnames(y)<-rep(g,2)
  return(y)
}

#logL
#f ... parameter to optimize
#fN ... list of cages to analyze (phenotype counts)
#fcm ... cas9 fitness mode 
#fom ... ot fitness mode 
#ft ... alread cut rate 
#fmpath ... path to transition matrix
#fmyg ... genotypes
#fx ... germline cut rate
#fe ... embryo cut-rate 
#return ... logLikelihood of data in N 

logL<-function(f,fN,fcm,fom,ft,fmpath,fmyg,fx,fe){
  #neutral set up:
  fmfcs<-1
  fmfch<-1
  fvcs<-1
  fvch<-1
  fmfos<-1
  fmfoh<-1
  fvos<-1
  fvoh<-1
  fne<-1
  fne<-f["ne"]
  
  fm<-init_m(fmpath,fmyg)
  fmmod<-mod_transition_m(fm,fx,fe,fmyg) #modified transition matrix
  
  #CAS9 
  if(fcm[1]>0) fmfcs<-f["mfcs"]
  if(fcm[1]==3) fmfch<-f["mfch"]
  if(fcm[2]>0) fmfcs<-f["mfcs"]
  if(fcm[2]==3) fmfch<-f["mfch"]
  if(fcm[3]>0) fvcs<-f["vcs"]
  if(fcm[3]==3) fvch<-f["vch"]
  
  #OT
  if(fom[1]>0) fmfos<-f["mfos"]
  if(fom[1]==3) fmfoh<-f["mfoh"]
  if(fom[2]>0) fmfos<-f["mfos"]
  if(fom[2]==3) fmfoh<-f["mfoh"]
  if(fom[3]>0) fvos<-f["vos"]
  if(fom[3]==3) fvoh<-f["voh"]
  
  fparams<-c(fne,fmfcs,fmfch,fvcs,fvch,fmfos,fmfoh,fvos,fvoh)
  names(fparams)<-c("ne","mfcs","mfch","vcs","vch","mfos","mfoh","vos","voh")
  
  l<-0 #start logL
  
  fn_sets<-length(fN)
  for(j in c(1:fn_sets)){
    fNj<-fN[[j]] #store counts
    fNfreq<-normalize_counts(fNj) #normalize counts
    fpar_geno<-start_genotype(x = fNfreq[1,],fg = fmyg,t = ft[j]) #initialize
    
    for(i in c(2:nrow(fNj))){
      fpro_pheno<-unlist(fNfreq[i,]) #obs phenotypes
      fpro_geno_exp<-propagate_genotypes(p=fpar_geno,m=fmmod,g=fmyg,f=fparams,fcas9 = fcm,fot = fom) #exp genotypes
      fpro_pheno_exp<-geno_to_pheno(fpro_geno_exp) #exp phenotypes given exp genotypes
      l<-l+rho(fpro_pheno,fpro_pheno_exp,fne) #add l
      fpar_geno<-pheno_to_geno(p=fpro_pheno,g=fpro_geno_exp) #obs genotypes i-1
    }
  }
  return(l)
}

#flnL ... log likelihood
#p ... number of params
#n ... number of data points
calculate_AICc<-function(flnL,p,n){
  aic=2*p-2*flnL+c(2*p*p+2*p)/c(n-p-1)
  return(aic)
}

#fml ... ML result 
#fdist ... min. dist for CI search
#falpha ... type 1 error
#fN ... effective population size 
#fcas9_modes ... cas9 fitness modes 
#fot_modes .... off-target fitness modes 
#g .... genotypes
#t ... already cut rate
#m ... matrix path 
#x ... germline cut rate 
#e ... embryo cut rate
calculate_CI<-function(fml,fdist,falpha,fN,fcas9_modes,fot_modes,g,t,m,x,e){
  fd<-fdist
  fsig<-qchisq(1-falpha,1) #calculate Chi^2 stat
  fLmax<-fml$value 
  y<-data.frame(par=fml$par)
  y$low_ci<-c(-1)
  y$up_ci<-c(-1)
  
  for(i in c(1:length(fml$par))){
    fdist<-ifelse(i>1,fd,1) #ne: other step limit
    fpara_dynamic<-fml$par #init
    fpara_test<-fpara_dynamic[i]
    fstep<-fpara_test/2
    frun<-T
    while(frun){
      fpara_test<-max(fdist,fpara_test-fstep)
      if(fpara_test==fdist) frun<-F #avoid endless loop
      fpara_dynamic[i]<-fpara_test #shift parameter 
      fl<-logL(f=fpara_dynamic,fN = fN,fcm = fcas9_modes,fom = fot_modes,ft=t,fmpath = m,fmyg = g,fx = x,fe = e) 
      fstat<-2*(fLmax-fl) #calculate Chi^2
      if(fstat>fsig){ #if significantly different
        if(2*fstep<fdist) {
          frun<-F
          y$low_ci[i]<-fpara_test
        } else {
          fpara_test<-fpara_test+fstep #restore parameter
          fstep<-fstep/2 #lower step size
        }
      }
    }
    
    fpara_dynamic<-fml$par #init upstream 
    fpara_test<-fpara_dynamic[i]
    fstep<-fpara_test/2
    frun<-T
    
    while(frun){
      fpara_test<-fpara_test+fstep #shift parameter
      if((fpara_test>=50000 & i==1) | (fpara_test>=5 & i>1)) run<-F #avoid endless loop 
      fpara_dynamic[i]<-fpara_test
      fl<-logL(f=fpara_dynamic,fN = fN,fcm = fcas9_modes,fom = fot_modes,ft=t,fmpath = m,fmyg = g,fx = x,fe = e) 
      fstat<-2*(fLmax-fl) #calculate Chi^2
      if(fstat>fsig){ #if significantly different
        if(2*fstep<fdist) {
          frun<-F
          y$up_ci[i]<-fpara_test
        } else {
          fpara_test<-fpara_test-fstep #restore parameter
          fstep<-fstep/2 #lower step size
        }
      }
    }
  }
  return(y)
}

#fcm ... cas9 fitness mode
#fom ... off-target fitness mode 
#fcm_limit ... upper cas9 param limit
#fom_limit ... upper ott-target param limit 
determine_para_range<-function(fcm,fom,fcm_limit=NULL,fom_limit=NULL){ #not all cross-checks of what is allowed 
  #full pstart
  fnames_param<-c("ne","mfcs","mfch","vcs","vch","mfos","mfoh","vos","voh")
  fpfull<-c(100,rep(1,8))
  names(fpfull)<-fnames_param
  fpfull_low<-c(25,rep(0.01,8))
  names(fpfull_low)<-fnames_param
  fpfull_up<-c(50000,rep(2,8))
  names(fpfull_up)<-fnames_param
  if(!is.null(fcm_limit)) fpfull_up[c("mfcs","vcs")]<-fcm_limit
  if(!is.null(fom_limit)) fpfull_up[c("mfos","vos")]<-fom_limit
  n<-c("ne")
  #CAS9 
  if(fcm[1]>0) n<-c(n,"mfcs")
  if(fcm[1]==3) n<-c(n,"mfch")
  if(fcm[2]>0) n<-c(n,"mfcs")
  if(fcm[2]==3) n<-c(n,"mfch")
  if(fcm[3]>0) n<-c(n,"vcs")
  if(fcm[3]==3) n<-c(n,"vch")
  
  #OT
  if(fom[1]>0) n<-c(n,"mfos")
  if(fom[1]==3) n<-c(n,"mfoh")
  if(fom[2]>0) n<-c(n,"mfos")
  if(fom[2]==3) n<-c(n,"mfoh")
  if(fom[3]>0) n<-c(n,"vos")
  if(fom[3]==3) n<-c(n,"voh")
  
  n<-n[!duplicated(n)]
  fpstart<-fpfull[n]
  fplow<-fpfull_low[n]
  fpup<-fpfull_up[n]
  return(list(fpstart,fpup,fplow))
}