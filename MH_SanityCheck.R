#### libraries loaded ####
library(ape)
library(phylobase)
library(mvtnorm)
library(MASS)

source("./Downloads/MH_Func.R")

#### hyperParam  ####
args <- commandArgs(TRUE)
print(args)
c_rate=c_shape=2
s_rate=s_shape=1
c_vec<-as.integer(c(seq(0.3,1,0.1),seq(1.5,3,0.5))*10)/10
sigma_vec<-seq(0.5,2,0.5)


obsData_idx<-as.numeric(args[1])
true_c<-as.integer(as.numeric(args[2])*10)/10
true_sigma2<-as.integer(as.numeric(args[3])*10)/10
chain_idx<-as.numeric(args[4])
c_idx<-which(c_vec==true_c)
s_idx<-which(sigma_vec==true_sigma2)

N=100
d=10

setwd('~/Research_Veera/PYDT/MCMC_result/DDT_MH')

#### Obs Data  ####

set.seed(202014063+ceiling(obsData_idx*100000)+ceiling(chain_idx*10000)+ceiling(c_idx*100)+ceiling(s_idx*100))

obsPhylo4d_lt<-readRDS(paste0("~/Downloads/obsData/obs_phylo4d_lt",obsData_idx,".RDS"))
#obsPhylo4d_lt<-readRDS(paste0("~/Research_Veera/PYDT/ABC_simData/DDT/N50d10/obsData/obs_phylo4d_lt1.RDS"))
obsDf<-obsPhylo4d_lt[[(c_idx-1)*length(sigma_vec)+s_idx]]@data[1:N,1:d]
rownames(obsDf)<-paste0("N",1:N)
if(as.integer(obsPhylo4d_lt[[(c_idx-1)*length(sigma_vec)+s_idx]]@data[1,"c"]*10)/10!=true_c | 
   as.integer(obsPhylo4d_lt[[(c_idx-1)*length(sigma_vec)+s_idx]]@data[1,"sigma2"]*10)/10!=true_sigma2){
  stop("different true param with df")
}

#### Init ####
options(warn=-1)

c_old<-true_c
s2_old<-true_sigma2
obs_hclu_phylo4d<-phylo4d(extractTree(obsPhylo4d_lt[[(c_idx-1)*length(sigma_vec)+s_idx]]),tip.data=obsDf)


#### MH Alg ####

#MCMC_num<-30000
#MCMC_num<-10000
MCMC_num<-5000
#MCMC_num<-3
print(paste0("total MCMC number: ",MCMC_num))
llh_chk<-data.frame(idx=numeric(0),
                    Accept=numeric(0),
                    llh_Prop=numeric(0),
                    q_ToOld=numeric(0),
                    llh_old=numeric(0),
                    q_ToNew=numeric(0),
                    c=numeric(0),
                    sigma2=numeric(0))
tree_lt<-list()
MRCA_mat_lt<-list()
trStruct_lt<-list()
options(warn=-1)



system.time({
  for(i in 1:MCMC_num){
    print(i)
    if(i==1){
      old_phylo4d<-obs_hclu_phylo4d
    }
    old_phylo4<-extractTree(old_phylo4d)
    rnd_Detach<-Random_Detach(tr_phylo4 = old_phylo4)
    prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                 old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                 c_ = c_old, theta_=0, alpha_=0)
    while(nTips(prop_newTree$new_phylo4)!=N){
      prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                   old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                   c_ = c_old, theta_=0, alpha_=0)
    }
    q_prob<-prop_log_prob(orgTree_phylo4 = old_phylo4, rmnTree = rnd_Detach$rmnTr,
                          old_detach_PaTime = rnd_Detach$paDivT, old_detach_PaLab = rnd_Detach$paNodeLab, old_detach_divLab = rnd_Detach$divLab,
                          new_divTime = prop_newTree$new_divT, new_AttachRoot = prop_newTree$attach_root, new_AttachTo = prop_newTree$attach_to,
                          c_=c_old)
    
    new_phylo4d<-phylo4d(prop_newTree$new_phylo4,tip.data=obsDf)

    if(i==1){
      oldTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = old_phylo4d)
      newTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = new_phylo4d)
    }else{
      oldTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = old_phylo4d,
                                trStruct_old = trStruct_lt[[(i-1)]] , tre_sigma_old = MRCA_mat_lt[[(i-1)]])
      newTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = new_phylo4d)
    }
    
    while(length(newTree_info)==1 | nTips(prop_newTree$new_phylo4)!=N){
      prop_newTree<-attach_subTree(subTree = rnd_Detach$subTr, rmnTree = rnd_Detach$rmnTr,
                                   old_divTime = rnd_Detach$divT, old_subTr_paLab = rnd_Detach$paNodeLab, 
                                   c_ = c_old, theta_=0, alpha_=0)
      q_prob<-prop_log_prob(orgTree_phylo4 = old_phylo4, rmnTree = rnd_Detach$rmnTr,
                            old_detach_PaTime = rnd_Detach$paDivT, old_detach_PaLab = rnd_Detach$paNodeLab, old_detach_divLab = rnd_Detach$divLab,
                            new_divTime = prop_newTree$new_divT, new_AttachRoot = prop_newTree$attach_root, new_AttachTo = prop_newTree$attach_to,
                            c_=c_old)
      
      new_phylo4d<-phylo4d(prop_newTree$new_phylo4,tip.data=obsDf)
      
      if(i==1){
        oldTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = old_phylo4d)
        newTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = new_phylo4d)
      }else{
        oldTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = old_phylo4d,
                                  trStruct_old = trStruct_lt[[(i-1)]] , tre_sigma_old = MRCA_mat_lt[[(i-1)]])
        newTree_info<-DDT_log_llh(c_=c_old, cmn_sigma2 = s2_old, tr_phylo4d = new_phylo4d)
      }
    }
    
    llh_chk[i,-2]<-c(i,newTree_info$res_log_llh,q_prob["q_RmnToOld"],
                     oldTree_info$res_log_llh,q_prob["q_RmnToNew"],
                     c_old, s2_old)
    u<-runif(1)
    alpha_Accp<-llh_chk[i,"llh_Prop"]+llh_chk[i,"q_ToOld"]-llh_chk[i,"llh_old"]-llh_chk[i,"q_ToNew"]
    if(u<exp(alpha_Accp)){
      llh_chk[i,2]<-1
      tree_lt[[i]]<-new_phylo4d
      MRCA_mat_lt[[i]]<-newTree_info$tre_sigma
      trStruct_lt[[i]]<-newTree_info$trStruct
      old_phylo4d<-new_phylo4d
      
    }else{
      llh_chk[i,2]<-0
      tree_lt[[i]]<-old_phylo4d
      MRCA_mat_lt[[i]]<-oldTree_info$tre_sigma
      trStruct_lt[[i]]<-oldTree_info$trStruct
    }
    
    c_old<-postSamp_c(shape_prior = c_shape, rate_prior = c_rate, trStruct_old = trStruct_lt[[i]])
    s2_old<-1/postSamp_s(shape_prior = s_shape, rate_prior = s_rate, MRCA_mat_old = MRCA_mat_lt[[i]])
    
  }
})


saveRDS(llh_chk,paste0("./MCMC_Compare/N100d10_Replicate/SanityCheck/llh_lt",obsData_idx,"_c",true_c,"_s",true_sigma2,"_ch",chain_idx,".RDS"))
saveRDS(tree_lt,paste0("./MCMC_Compare/N100d10_Replicate/SanityCheck/tree_lt",obsData_idx,"_c",true_c,"_s",true_sigma2,"_ch",chain_idx,".RDS"))
options(warn=0)
