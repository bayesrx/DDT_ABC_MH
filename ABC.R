setwd('~/Research_Veera/PYDT/ABC_simData/DDT/N100d10/')
library(abc)
library(ape)
library(phylobase)

#### Load sim data ####
N=100
d=10
f_names<-sort(list.files("./")[grepl("sim_sumstat",list.files())])

for(i in 1:length(f_names)){
  print(i)
  f_n<-f_names[i]
  if(i==1){
    sim_sumstat<-readRDS(paste0("./",f_n))
  }else{
    ele<-readRDS(paste0("./",f_n))
    sim_sumstat<-rbind(sim_sumstat,ele)
  }
}
rm(ele)

na_idx<-apply(sim_sumstat,1,function(x) any(is.na(x)))
if(any(na_idx)){
  sim_sumstat<-sim_sumstat[!na_idx,]
}

sim_param<-data.frame(sim_sumstat[,c("idx","c_hazard","sigma_sq")])
sim_sumstat<-data.frame(sim_sumstat[,c(paste0("dist_org",1:N),paste0("dist_rel",1:99),paste0("hclu_tip_BrLen",c(10,25,50,75,90)))])
sim_sumstat$org2_mean<-apply(sim_sumstat[,paste0("dist_org",1:N)],1,function(x) sum(x^2)/(d*N))


data_lt_vec<-101:130


#### Obs Data ####
c_vec<-as.integer(c(seq(0.3,1,0.1),seq(1.5,3,0.5))*10)/10
sigma_vec<-seq(0.5,2,0.5)
options(warn=-1)

obs_sumstat<-data.frame(matrix(NA,nrow=length(c_vec)*length(sigma_vec)*length(data_lt_vec), ncol=(3+N+99+5+1)))
colnames(obs_sumstat)<-c("data","c_hazard","sigma_sq",
                         paste0("dist_org",1:N),
                         paste0("dist_rel",1:99),
                         paste0("hclu_tip_BrLen",c(10,25,50,75,90)),
                         "org2_mean")
for(data_idx in 1:length(data_lt_vec)){
  obsData_idx<-data_lt_vec[data_idx]
  print(obsData_idx)
  obsPhylo4d_lt<-readRDS(paste0("./obsData/obs_phylo4d_lt",obsData_idx,".RDS"))
  for(c_idx in 1:length(c_vec)){
    for(s_idx in 1:length(sigma_vec)){
      idx<-(data_idx-1)*length(c_vec)*length(sigma_vec)+(c_idx-1)*length(sigma_vec)+s_idx
      obsDf<-obsPhylo4d_lt[[(c_idx-1)*length(sigma_vec)+s_idx]]@data[1:N,]
      rownames(obsDf)<-paste0("N",1:N)
      obs_sumstat[idx,c("data","c_hazard","sigma_sq")]<-c(obsData_idx,obsDf[1,c("c","sigma2")])
      tip_loc<-obsDf[,1:d]
      tip_dist<-dist(tip_loc)
      dist_org<-sort(apply(tip_loc,1,function(x) sqrt(sum(x^2))))
      dist_rel_pct<-quantile(as.vector(tip_dist),seq(0.01,0.99,0.01))
      hclu<-as.phylo(hclust(tip_dist))
      hclu_tip_BrLen<-quantile(edgeLength(phylo4(hclu))[getEdge(phylo4(hclu),1:N,"descendant")],c(0.1,0.25,0.5,0.75,0.9))
      obs_sumstat[idx,c(paste0("dist_org",1:N),paste0("dist_rel",1:99),paste0("hclu_tip_BrLen",c(10,25,50,75,90)),"org2_mean")]<-
        c(dist_org,dist_rel_pct,hclu_tip_BrLen,sum(dist_org^2)/(d*N))
    }
  }
}




#### ABC Inference ####
set.seed(1234)
abc_s_adj<-list()
abc_c_adj<-list()
abc_s_w<-list()
abc_c_w<-list()


NSS<-3000
totalSimNum<-nrow(sim_sumstat)
frac<-totalSimNum/nrow(sim_sumstat)
sim_idx<-sort(sample(1:nrow(sim_sumstat),nrow(sim_sumstat)*frac))
sim_sumstat_frac<-sim_sumstat[sim_idx,]
sim_param_frac<-sim_param[sim_idx,]

for(data_idx in 1:length(data_lt_vec)){
  obsData_idx<-data_lt_vec[data_idx]
  print(obsData_idx)
  abc_c_summary<-abc_s_summary<-matrix(NA,ncol=7,nrow=8)
  colnames(abc_s_summary)<-colnames(abc_c_summary)<-c("data","c","s2","mean","2.5P","97.5P","median")
  lt_idx<-1
  count_idx<-1
  for(c_idx in 1:length(c_vec)){
    c_hazard=c_vec[c_idx]
    if(!(c_hazard %in% c(0.3,0.5,0.7,1))){
      next
    }
    #print(c_hazard)
    for(s_idx in 1:length(sigma_vec)){
      sigma2<-sigma_vec[s_idx]
      if(!(sigma2 %in% c(0.5,1))){
        next
      }
      #print(sigma2)
      idx<-(data_idx-1)*length(c_vec)*length(sigma_vec)+(c_idx-1)*length(sigma_vec)+s_idx
      abc_s<-abc(target=as.matrix(obs_sumstat[idx,"org2_mean"]),
                 sumstat=as.matrix(sim_sumstat_frac[,"org2_mean"]),
                 param=sim_param_frac[,c("sigma_sq")],
                 tol=NSS/totalSimNum,transf = rep("none",1),method ="loclinear")
      
      abc_c<-abc(target=as.matrix(obs_sumstat[idx,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),
                                                    paste0("dist_rel",c(10,25,50,75,90)))]),
                 sumstat=sim_sumstat_frac[,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),
                                        paste0("dist_rel",c(10,25,50,75,90)))],
                 param=sim_param_frac$c_hazard,
                 tol=NSS/totalSimNum,transf = rep("none",1),method ="loclinear")
      
      abc_s_adj[[lt_idx]]<-abc_s$adj.values
      abc_c_adj[[lt_idx]]<-abc_c$adj.values
      abc_s_w[[lt_idx]]<-abc_s$weights
      abc_c_w[[lt_idx]]<-abc_c$weights
      names(abc_s_adj)[lt_idx]<-names(abc_c_adj)[lt_idx]<-names(abc_s_w)[lt_idx]<-names(abc_c_w)[lt_idx]<-paste0("c",c_hazard,"s",sigma2)
      tmp_s<-summary(abc_s,print=F)
      tmp_c<-summary(abc_c,print=F)
      abc_s_summary[count_idx,]<-c(obsData_idx,c_hazard,sigma2,tmp_s[c(4,2,6,5),])
      abc_c_summary[count_idx,]<-c(obsData_idx,c_hazard,sigma2,tmp_c[c(4,2,6,5),])
      lt_idx<-lt_idx+1
      count_idx<-count_idx+1
    }
  }
  saveRDS(list(abc_s_adj=abc_s_adj,abc_c_adj=abc_c_adj,
               abc_s_w=abc_s_w,abc_c_w=abc_c_w,
               abc_c_summary=abc_c_summary,abc_s_summary=abc_s_summary),paste0("./ABC_results/wtd_frac_lt",obsData_idx,".RDS"))
}

#### time control ####
set.seed(12345)
abc_s_adj<-list()
abc_c_adj<-list()
abc_s_w<-list()
abc_c_w<-list()
totalSimNum<-round(3.2872 / 1.8952 * 10000)
frac<-totalSimNum/nrow(sim_sumstat)
sim_idx<-sort(sample(1:nrow(sim_sumstat),nrow(sim_sumstat)*frac))
sim_sumstat_frac<-sim_sumstat[sim_idx,]
sim_param_frac<-sim_param[sim_idx,]
threshold_<-0.05

for(data_idx in 1:length(data_lt_vec)){
  obsData_idx<-data_lt_vec[data_idx]
  print(obsData_idx)
  lt_idx<-1
  abc_c_summary<-abc_s_summary<-matrix(NA,ncol=7,nrow=8)
  colnames(abc_s_summary)<-colnames(abc_c_summary)<-c("data","c","s2","mean","2.5P","97.5P","median")
  for(c_idx in 1:length(c_vec)){
    c_hazard=c_vec[c_idx]
    if(!(c_hazard %in% c(0.3,0.5,0.7,1))){
      next
    }
    print(c_hazard)
    for(s_idx in 1:length(sigma_vec)){
      sigma2<-sigma_vec[s_idx]
      if(!(sigma2 %in% c(0.5,1))){
        next
      }
      idx<-(data_idx-1)*length(c_vec)*length(sigma_vec)+(c_idx-1)*length(sigma_vec)+s_idx
      abc_s<-abc(target=as.matrix(obs_sumstat[idx,"org2_mean"]),
                 sumstat=as.matrix(sim_sumstat_frac[,"org2_mean"]),
                 param=sim_param_frac[,c("sigma_sq")],
                 tol=threshold_,transf = rep("none",1),method ="loclinear")
      
      abc_c<-abc(target=as.matrix(obs_sumstat[idx,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),
                                                    paste0("dist_rel",c(10,25,50,75,90)))]),
                 sumstat=sim_sumstat_frac[,c(paste0("hclu_tip_BrLen",c(10,25,50,75,90)),
                                             paste0("dist_rel",c(10,25,50,75,90)))],
                 param=sim_param_frac$c_hazard,
                 tol=threshold_,transf = rep("none",1),method ="loclinear")
      
      abc_s_adj[[lt_idx]]<-abc_s$adj.values
      abc_c_adj[[lt_idx]]<-abc_c$adj.values
      abc_s_w[[lt_idx]]<-abc_s$weights
      abc_c_w[[lt_idx]]<-abc_c$weights
      names(abc_s_adj)[lt_idx]<-names(abc_s_adj)[lt_idx]<-names(abc_s_w)[lt_idx]<-names(abc_c_w)[lt_idx]<-paste0("c",c_hazard,"s",sigma2)
      tmp_s<-summary(abc_s,print=F)
      tmp_c<-summary(abc_c,print=F)
      abc_s_summary[lt_idx,]<-c(obsData_idx,c_hazard,sigma2,tmp_s[c(4,2,6,5),])
      abc_c_summary[lt_idx,]<-c(obsData_idx,c_hazard,sigma2,tmp_c[c(4,2,6,5),])
      lt_idx<-lt_idx+1
    }
  }
  saveRDS(list(abc_s_adj=abc_s_adj,abc_c_adj=abc_c_adj,
               abc_s_w=abc_s_w,abc_c_w=abc_c_w,
               abc_c_summary=abc_c_summary,abc_s_summary=abc_s_summary),
          paste0("./ABC_results/wtd_tCtrl_",threshold_*100,"Per_lt",obsData_idx,".RDS"))
}


