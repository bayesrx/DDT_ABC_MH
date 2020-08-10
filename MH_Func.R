#### libraries loaded ####
library(ape)
library(phylobase)
library(mvtnorm)
library(MASS)

#### User Defined Function  ####
is.invertible <- function(m){
  class(try(solve(m),silent=T))=="matrix"
} 

Marginal_loc_log_llh<-function(tr_Phylo4d,cmn_Sigma2,tre_sigma_old){
  loc_<-tr_Phylo4d@data[1:N,1:d]
  rownames(loc_)<-tipLabels(tr_Phylo4d)
  loc_<-loc_[order(as.numeric(substring(rownames(loc_),2))),]
  if(missing(tre_sigma_old)){
    tre_sigma<-matrix(NA,nrow=N,ncol=N)
    N_grid<-expand.grid(1:N,1:N)
    N_grid<-N_grid[N_grid$Var1>=N_grid$Var2,]
    tre_sigma[lower.tri(tre_sigma,diag = T)]<-apply(N_grid,1,function(x) nodeHeight(tr_Phylo4d,MRCA(tr_Phylo4d,paste0("N",x[1]),paste0("N",x[2])),"root"))
    tre_sigma[upper.tri(tre_sigma)]<-t(tre_sigma)[upper.tri(tre_sigma)]
  }else{
    tre_sigma=tre_sigma_old
  }
  if(!is.invertible(tre_sigma)){
    return(list(singularMat=TRUE))
  }else{
    return(list(llh=sum(mvtnorm::dmvnorm(t(loc_),mean=rep(0,N),sigma=cmn_Sigma2*tre_sigma,log=T)),MRCA_mat=tre_sigma))
  }
}

H_n<-function(n_){
  if(n_==0){
    return(0)
  }else{
    return(sum(1/1:n_))
  }
}

J_n<-function(l_v,r_v){
  return(H_n(n_=(r_v+l_v-1))-H_n(n_=(l_v-1))-H_n(n_=(r_v-1)))
}

divTime_DDT_log_llh<-function(c_,l_v,r_v,t_v){
  J=J_n(l_v=l_v,r_v=r_v)
  return(log(c_)+(c_*J-1)*log(1-t_v))
}

treeStruct_DDT_log_llh<-function(l_v,r_v){
  return(lfactorial(l_v-1)+ lfactorial(r_v-1) - lfactorial(r_v+l_v-1))
}

DDT_log_llh<-function(c_,cmn_sigma2,tr_phylo4d,trStruct_old,tre_sigma_old){
  if(missing(trStruct_old)){
    tr_phylo4<-extractTree(tr_phylo4d)
    all_divT<-c(nodeHeight(tr_phylo4,nodeLabels(tr_phylo4),from="root"),rep(1,N))
    names(all_divT)[(N+1):(2*N)]<-1:N
    from_NumPt<-c(sapply(descendants(tr_phylo4,nodeLabels(tr_phylo4),"tip"),length),rep(1,N))
    names(from_NumPt)[(N+1):(2*N)]<-1:N
    trStruct_df<-data.frame(node=names(all_divT),
                            divT=all_divT,
                            m_v=from_NumPt)
    all_edge<-getEdge(tr_phylo4)
    allEdge_dashIdx<-unlist(gregexpr(pattern="-",all_edge))
    trStruct<-data.frame(from_=as.numeric(substring(all_edge,1,allEdge_dashIdx-1)),
                         to_=as.numeric(substring(all_edge,allEdge_dashIdx+1)))
    trStruct<-trStruct[trStruct$from_!=0,]
    trStruct<-trStruct[trStruct$to_ >N,]
    trStruct<-merge(x=trStruct,y=trStruct_df[,c("node","divT")],by.x="to_",by.y="node",all.x=T)
    trStruct<-merge(x=trStruct,y=trStruct_df[,c("node","m_v")],by.x="to_",by.y="node",all.x=T)
    trStruct_dup<-merge(x=trStruct,y=trStruct[,c("from_","m_v")],by.x="to_",by.y="from_",all.x=T,all.y=F)
    trStruct_dup<-trStruct_dup[!duplicated(trStruct_dup[,1:3]),]
    colnames(trStruct_dup)[4:5]<-c("m_v","l_v")
    trStruct_dup$l_v[trStruct_dup$m_v==2]=1
    trStruct_dup$r_v<-trStruct_dup$m_v-trStruct_dup$l_v
    trStruct_dup$trStruct_llh<-apply(trStruct_dup,1,function(x) treeStruct_DDT_log_llh(l_v=x[5],r_v=x[6]))
    trStruct_dup$divT_llh<-apply(trStruct_dup,1,function(x) divTime_DDT_log_llh(c_=c_,l_v=x[5],r_v=x[6],t_v=x[3]))
  }else{
    trStruct_dup=trStruct_old
  }
  location_llh<-Marginal_loc_log_llh(tr_Phylo4d =  tr_phylo4d,cmn_Sigma2 =  cmn_sigma2, tre_sigma_old )
  if(length(location_llh)==1){
    return(location_llh)
  }else{
    res_log_llh<-sum(trStruct_dup$trStruct_llh,trStruct_dup$divT_llh) +location_llh$llh
    return(list(res_log_llh=res_log_llh, trStruct=trStruct_dup, tre_sigma=location_llh$MRCA_mat))
  }
}

a_t<-function(t,c){
  return(c/(1-t))
}

A_t<-function(t,c){
  return(-c*log(1-t))
}

A_inv<-function(A,c){
  return(1-exp(-A/c))
}

div_time<-function(m_v,t_u,c_,theta_,alpha_){
  u=runif(1)
  input=A_t(t_u, c=c_)-exp(lgamma(m_v+1+theta_)-lgamma(m_v-alpha_))*log(1-u)
  return(A_inv(input, c=c_))
}

mean_divTime<-function(a, t_u){
  return((1-pbeta(t_u,2,a))/((a+1)*(1-t_u)^a))
}

add_tip<-function(tr_old,div_t,to_time,to_node,label){
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=label,
            edge.length=1-div_t,
            Nnode=1)
  class(tip)<-"phylo"
  new_tr<-bind.tree(tr_old,tip,where=to_node,position=to_time-div_t)
  return(new_tr)
}

add_root<-function(tr,root_edge,rootLabel,tipLabel){
  tip<-list(edge=matrix(c(2,1),1,2),
            node.label=rootLabel,
            tip.label=tipLabel,
            edge.length=root_edge,
            Nnode=1)
  class(tip)<-"phylo"
  rooted_Ape<-bind.tree(tip,tr,where=1)
  return(rooted_Ape)
}

Random_Detach<-function(tr_phylo4){
  root_name<-names(rootNode(tr_phylo4))
  detachNode<-sample( setdiff(paste0("N",(1:(nNodes(tr_phylo4)+nTips(tr_phylo4)))),
                              c(root_name,names(descendants(tr_phylo4,root_name,"children")))),
                      1,F)
  paNode<-ancestor(tr_phylo4,detachNode)
  paEdge<-edgeLength(tr_phylo4)[getEdge(tr_phylo4,detachNode,"descendant")]
  paDivTime<-nodeHeight(tr_phylo4,paNode,"root")
  divTime<-nodeHeight(tr_phylo4,detachNode,"root")
  subTree_tip<-names(descendants(tr_phylo4,detachNode,"tips"))
  #if(names(paNode)==paste0("N",N+2)){}
  if(length(subTree_tip)==1){
    subTree<-detachNode
  }else{
    subTree<-as(subset(tr_phylo4,node.subtree=detachNode),"phylo")  
  }
  # remaining tree contains the singleton from thr original root 
  rmnTree_tip<-setdiff(paste0("N",1:N),subTree_tip)
  if(length(rmnTree_tip)==1){
    rmnTree<-list(edge=matrix(c(2,1),1,2),
                  node.label=paste0("N",N+1),
                  tip.label=rmnTree_tip,
                  edge.length=1,
                  Nnode=1)
    class(rmnTree)<-"phylo"
  }else{
    rmnTree<-add_root(tr=as(subset(tr_phylo4,tips.include=rmnTree_tip),"phylo"),
                      root_edge = nodeHeight(tr_phylo4,names(MRCA(tr_phylo4,rmnTree_tip)),from="root"),
                      rootLabel = paste0("N",N+1),
                      tipLabel = names(rootNode(tr_phylo4)))
  }
  return(list(subTr=subTree,rmnTr=rmnTree,paNodeLab=names(paNode),paDivT=paDivTime,divT=divTime, divLab=detachNode))
}

ReAttach_Pt<-function(rmnTree,divTime,c_,theta_,alpha_){
  rmnTree_tip<-Ntip(rmnTree)
  new_tip<-list()
  if(rmnTree_tip==1){
    t<-div_time(1,0,c_=true_c ,theta_ = theta, alpha_= alpha )
    return(list(divT=t,divRoot=rmnTree$node.label,divTo=rmnTree$tip.label,distToDiv=1-t))
  }else{
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    Cur_n<-rmnTree_tip
    root_time=0
    root_<-names(rootNode(rmnTree_phylo4))
    root_to<-names(descendants(rmnTree_phylo4, root_, "child"))
    edge_len<-edgeLength(rmnTree_phylo4)[getEdge(rmnTree_phylo4, root_to)]
    dist_root_to<-root_time+edge_len
    while(length(new_tip)==0){
      t<-div_time(Cur_n,root_time, c_ ,theta_, alpha_)
      if(t<dist_root_to){
        return(list(divT=t,divRoot=root_,divTo=root_to,distToDiv=dist_root_to-t))
      }
      root_div<-names(descendants(rmnTree_phylo4, root_to, type="child"))
      K=length(root_div)
      n_k<-unlist(lapply(descendants(rmnTree_phylo4, root_div, type="tip"),length))
      names(n_k)<-root_div
      prob_vec<-c(n_k-alpha_,theta_+alpha_*K)/(sum(n_k)+theta_)
      path_idx_sel<-which.max(rmultinom(1,1,prob_vec))
      if(path_idx_sel==length(prob_vec)){
        return(list(divT=dist_root_to,divRoot=root_,divTo=root_to,distToDiv=dist_root_to-t))
      }else{
        selected_node<-root_div[path_idx_sel]
        m_v<-n_k[path_idx_sel]
        if(m_v==1){
          t<-div_time(1,dist_root_to, c_,theta_, alpha_)
          return(list(divT=t,divRoot=root_to,divTo=selected_node,distToDiv=1-t))
        }else{
          root_time=dist_root_to
          Cur_n<-m_v
          root_<-root_to
          root_to<-selected_node
          edge_len<-edgeLength(rmnTree_phylo4)[getEdge(rmnTree_phylo4, root_to)]
          dist_root_to<-root_time+edge_len
        }
      }
    }
  }
}

attach_subTree<-function(subTree, rmnTree, old_divTime, old_subTr_paLab, c_, theta_, alpha_){
  if(is.character(subTree)){
    newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    while(newDiv_pt$divT>=old_divTime){
      newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    }
    new_phylo4<-phylo4(add_tip(rmnTree,
                               div_t = newDiv_pt$divT, 
                               to_time = nodeHeight(rmnTree_phylo4,newDiv_pt$divTo,"root"), 
                               to_node=getNode(rmnTree_phylo4,newDiv_pt$divTo),
                               label = subTree))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=newDiv_pt$divRoot,attach_to=newDiv_pt$divTo,new_divT=newDiv_pt$divT))
  }else if(Ntip(rmnTree)==1){
    t<-div_time(1,0,c_=c_ ,theta_ = theta_, alpha_= alpha_ )
    while(t>=old_divTime){
      t<-div_time(1,0,c_=c_ ,theta_ = theta_, alpha_= alpha_ )
    }
    new_subTree<-add_root(tr=subTree,
                          root_edge = old_divTime - t ,
                          rootLabel = old_subTr_paLab,
                          tipLabel = names(rootNode(phylo4(subTree))))
    
    new_phylo4<-phylo4(bind.tree(rmnTree,
                                 new_subTree,
                                 where=1, 
                                 position = 1-t))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=rmnTree$node.label,attach_to=rmnTree$tip.label,new_divT=t))
  }else{
    newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    rmnTree_phylo4<-phylo4(rmnTree,check.node.labels="keep")
    while(newDiv_pt$divT>=old_divTime){
      newDiv_pt<-ReAttach_Pt(rmnTree=rmnTree, c_=c_, theta_=theta_, alpha_=alpha_)
    }
    new_subTree<-add_root(tr=subTree,
                          root_edge = old_divTime - newDiv_pt$divT ,
                          rootLabel = old_subTr_paLab,
                          tipLabel = names(rootNode(phylo4(subTree))))
    new_phylo4<-phylo4(bind.tree(rmnTree,
                                 new_subTree,
                                 where=getNode(rmnTree_phylo4,newDiv_pt$divTo), 
                                 position = newDiv_pt$distToDiv))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))]<-old_subTr_paLab
    return(list(new_phylo4=new_phylo4, attach_root=newDiv_pt$divRoot,attach_to=newDiv_pt$divTo,new_divT=newDiv_pt$divT))
  }
}

prop_log_prob<-function(orgTree_phylo4 ,rmnTree, old_detach_PaTime, old_detach_PaLab, old_detach_divLab, new_divTime, new_AttachRoot, new_AttachTo,c_){
  if(nTips(rmnTree)==1){
    q_RmnToNew=-A_t(new_divTime,c=c_)+log(a_t(new_divTime,c=c_))
    q_RmnToOld=-A_t(old_detach_PaTime,c=c_)+log(a_t(old_detach_PaTime,c=c_))
    return(c(q_RmnToNew=q_RmnToNew, q_RmnToOld=q_RmnToOld))
  }else{
    rmnTree_phylo4<-phylo4(rmnTree)
    rmnTree_root_<-names(rootNode(rmnTree_phylo4))
    
    old_DeTachRoot<-names(ancestor(orgTree_phylo4,old_detach_PaLab))
    old_DeTachTo<-setdiff(names(descendants(orgTree_phylo4,old_detach_PaLab,"children")),old_detach_divLab)
    old_pathRoot<-unique(c(rmnTree_root_,names(shortestPath(rmnTree_phylo4,old_DeTachRoot,rmnTree_root_)),old_DeTachRoot,old_DeTachTo))
    if(length(old_pathRoot)==2){
      m_v<-length(descendants(rmnTree_phylo4,old_pathRoot[-1],"tips"))
      names(m_v)<-old_pathRoot[-1]
      divT<-old_detach_PaTime
      names(divT)<-"old_detachPaTime"
    }else{
      m_v<-sapply(descendants(rmnTree_phylo4,old_pathRoot[-1],"tips"),length)
      names(m_v)<-old_pathRoot[-1]
      divT<-c(nodeHeight(rmnTree_phylo4,names(m_v)[-length(m_v)],"root"),old_detach_PaTime)
      names(divT)<-c(names(m_v)[-length(m_v)],"old_detachPaTime")
    }
    A_t_divT<-A_t(divT,c=c_)
    A_t_tu<-c(0,A_t_divT[-length(A_t_divT)])
    A_t_tv<-A_t_divT
    final_A_t<-sum((A_t_tu-A_t_tv)/m_v)
    frac<-sum(log(m_v[-1]/m_v[-length(m_v)]))
    final_a_t<-log(a_t(old_detach_PaTime,c=c_)/m_v[length(m_v)])
    q_RmnToOld=sum(final_A_t,frac,final_a_t)
    
    new_pathRoot<-unique(c(rmnTree_root_,names(shortestPath(rmnTree_phylo4,new_AttachRoot,rmnTree_root_)),new_AttachRoot,new_AttachTo))
    if(length(new_pathRoot)==2){
      m_v<-length(descendants(rmnTree_phylo4,new_pathRoot[-1],"tips"))
      names(m_v)<-new_pathRoot[-1]
      divT<-new_divTime
      names(divT)<-"new_divTime"
    }else{
      m_v<-sapply(descendants(rmnTree_phylo4,new_pathRoot[-1],"tips"),length)
      names(m_v)<-new_pathRoot[-1]
      divT<-c(nodeHeight(rmnTree_phylo4,names(m_v)[-length(m_v)],"root"),new_divTime)
      names(divT)<-c(names(m_v)[-length(m_v)],"new_divTime")
    }
    A_t_divT<-A_t(divT,c=c_)
    A_t_tu<-c(0,A_t_divT[-length(A_t_divT)])
    A_t_tv<-A_t_divT
    final_A_t<-sum((A_t_tu-A_t_tv)/m_v)
    frac<-sum(log(m_v[-1]/m_v[-length(m_v)]))
    final_a_t<-log(a_t(new_divTime,c=c_)/m_v[length(m_v)])
    q_RmnToNew=sum(final_A_t,frac,final_a_t)
    
    return(c(q_RmnToNew=q_RmnToNew, q_RmnToOld=q_RmnToOld))
  }
}

postSamp_c<-function(shape_prior,rate_prior,trStruct_old){
  new_shape=shape_prior+nrow(trStruct_old)
  new_rate=rate_prior + sum((sapply(trStruct_old$m_v-1,function(x) H_n(x)) - 
                               sapply(trStruct_old$l_v-1,function(x) H_n(x)) - 
                               sapply(trStruct_old$r_v-1,function(x) H_n(x))) * (-log(1- trStruct_old$divT)))
  return(rgamma(1,shape=new_shape,rate=new_rate))
}

postSamp_s<-function(shape_prior,rate_prior,MRCA_mat_old){
  new_shape = shape_prior + N*d/2
  inv_MRCA = solve(MRCA_mat_old)
  new_rate = rate_prior + sum(apply(obsDf,2,function(x) x %*% inv_MRCA %*% as.matrix(x,ncol=1)))/2
  return(rgamma(1,shape=new_shape,rate=new_rate))
}