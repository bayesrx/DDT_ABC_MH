library(ape)
library(phylobase)
library(ggplot2)
source("~/Downloads/ABC_Func.R")


args <- commandArgs(TRUE)
print(args)
args_seed<-as.numeric(args[1])



set.seed(args_seed)
#N=50
N=100
d=10

#sim_max<-5000
sim_max<-100
sim_phylo4<-list()
sim_sumstat<-data.frame(matrix(NA,nrow=sim_max, ncol=(5+N+99+5)))
colnames(sim_sumstat)<-c("idx","c_hazard","sigma_sq","alpha","theta",
                         paste0("dist_org",1:N),
                         paste0("dist_rel",1:99),
                         paste0("hclu_tip_BrLen",c(10,25,50,75,90)))



alpha_=0
theta_=0

options(warn=-1)


system.time({
  for(sim_idx in 1:sim_max){
    tryCatch({
      print(sim_idx)
      c_hazard<-rgamma(1,2,rate=2)
      sigma_sq<-1/rgamma(1,1,1)
      for(n in 2:N){
        if(n==2){
          t<-div_time(1,0,c_=c_hazard ,theta_ = theta_, alpha_=alpha_)
          tree_txt<-paste("((1:",1-t,",2:",1-t,"):",t,");",sep='')
          tr<-read.tree(text=tree_txt)
        }else{
          tr<-add_One_Obs(tr, c_=c_hazard, alpha_= alpha_, theta_=theta_)
        }
      }
      #plot(as(tr,"phylo4"))
      
      if(Ntip(tr)!=N) next
      tr_phylo4d<-gen_location_homeskd(tr, sigma2=diag(rep(sigma_sq,d)), dime=d)
      cIdx_df<-matrix(c(rep(c_hazard,N),rep(sim_idx,N)),ncol=2)
      colnames(cIdx_df)=c("c","idx")
      tr_phylo4d<-addData(tr_phylo4d,tip.data=cIdx_df)
      sim_phylo4[[sim_idx]]<-tr_phylo4d
      if(Ntip(tr)!=N){
        sim_sumstat[sim_idx,]<-c(sim_idx,c_hazard,sigma_sq,alpha_,theta_,rep(NA,ncol(sim_sumstat)-5))
      }else{
        tip_loc<-tr_phylo4d@data[1:N,1:d]
        tip_dist<-dist(tip_loc)
        dist_org<-sort(apply(tip_loc,1,function(x) sqrt(sum(x^2))))
        dist_rel_pct<-quantile(as.vector(tip_dist),seq(0.01,0.99,0.01))
        hclu<-as.phylo(hclust(tip_dist))
        #c_MME_hclu<-c_MME_cal(hclu)
        hclu_tip_BrLen<-quantile(edgeLength(phylo4(hclu))[getEdge(phylo4(hclu),1:N,"descendant")],c(0.1,0.25,0.5,0.75,0.9))
        #sim_sumstat[sim_idx,]<-c(sim_idx,c_hazard,sigma_sq,alpha_,theta_,c_MME_hclu,dist_org,dist_rel_pct,hclu_tip_brLen_quan)
        sim_sumstat[sim_idx,]<-c(sim_idx,c_hazard,sigma_sq,alpha_,theta_,dist_org,dist_rel_pct,hclu_tip_BrLen)
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
})

saveRDS(sim_phylo4,file=paste0("/home/yaots/Research_Veera/PYDT/ABC_simData/DDT/N100d10/sim_phylo4d_",args_seed,".RDS"))
saveRDS(sim_sumstat,file=paste0("/home/yaots/Research_Veera/PYDT/ABC_simData/DDT/N100d10/sim_sumstat_",args_seed,".RDS"))

