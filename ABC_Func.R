library(ape)
library(phylobase)
library(ggplot2)

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

add_multichotomous_tip<-function(tr_old,div_t,div_node,label){
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=label,
            edge.length=1-div_t,
            Nnode=1)
  class(tip)<-"phylo"
  new_tr<-bind.tree(tr_old,tip,where=div_node)
  return(new_tr)
}

add_One_Obs<-function(input_tr, c_, alpha_, theta_){
  n<-Ntip(input_tr)
  #print(paste("adding", n+1, "points"))
  tr<-input_tr
  tr$tip.label<-paste("N",as.character(1:n), sep="")
  tr$node.label<-paste("N",as.character((n+1):(n+Nnode(tr))), sep="")
  tr_phylo4<-phylo4(tr,check.node.labels="keep")
  dist_node<-nodeHeight(tr_phylo4, from="root")
  new_tip<-list()
  Cur_n<-n
  root_time=0
  root_<-rootNode(tr_phylo4)
  root_to<-descendants(tr_phylo4, root_, "child")
  edge_len<-edgeLength(tr_phylo4)[getEdge(tr_phylo4, names(root_to))]
  dist_root_to<-root_time+edge_len
  while(length(new_tip)==0){
    t<-div_time(Cur_n,root_time, c_ ,theta_, alpha_)
    #print(paste("root_: ", names(root_),"root_time: ", root_time, "root_to: ", names(root_to)))
    #print(paste("Cur_n: ", Cur_n,"t: ",t, "edge_len: ", edge_len, "dist(root, root_to): ", dist_root_to))
    if(t<dist_root_to){
      #print(paste("t<dist_node[as.numeric(root_to_label)]"))
      #print("")
      #print("")
      new_tr<-add_tip(tr,t,dist_root_to,
                      as.numeric(substr(names(root_to),2, nchar(names(root_to)))),
                      paste("N", as.character(n+1), sep=""))
      return(new_tr)
    }
    root_div<-descendants(tr_phylo4, root_to, type="child")
    #print(paste("root_div: ", names(root_div)))
    K=length(root_div)
    n_k<-unlist(lapply(descendants(tr_phylo4, names(root_div), type="tip"),length))
    prob_vec<-c(n_k-alpha_,theta_+alpha_*K)/(sum(n_k)+theta_)
    path_idx_sel<-which.max(rmultinom(1,1,prob_vec))
    if(path_idx_sel==length(prob_vec)){
      #print(paste("path_idx_sel==length(prob_vec);")) 
      #print(paste("root_to: ", names(root_to), "; dist(root, root_to): ", dist_root_to))
      #print("")
      #print("")
      new_tr<-add_multichotomous_tip(tr,dist_root_to,
                                     as.numeric(substr(names(root_to),2, nchar(names(root_to)))),
                                     paste("N",as.character(n+1),sep=""))
      return(new_tr)
    }else{
      selected_node<-root_div[path_idx_sel]
      #print(paste("selected_node: ", names(selected_node)))
      m_v<-n_k[path_idx_sel]
      if(m_v==1){
        t<-div_time(1,dist_root_to, c_,theta_, alpha_)
        #print(paste("m_v==1"))
        #print(paste("root_to: ", names(root_to), "; dist(root_to): ", dist_root_to, "t: ", t))
        #print("")
        #print("")
        new_tr<-add_tip(tr, t, 1, 
                        as.numeric(substr(names(selected_node),2, nchar(names(selected_node)))),
                        paste("N",as.character(n+1),sep=""))
        return(new_tr)
      }else{
        root_time=dist_root_to
        Cur_n<-m_v
        root_<-root_to
        root_to<-selected_node
        edge_len<-edgeLength(tr_phylo4)[getEdge(tr_phylo4, names(root_to))]
        dist_root_to<-root_time+edge_len
        #print(paste("root_time", root_time))
        #plot(cur_tr4,show.node.label = T)
      }
    }
  }
}

chol_PD<-function(mat){
  init<-tryCatch({
    base_chol<-chol(mat)
    noError<-T
  },error=function(e){
    print(e)
    mat<-diag(1/rgamma(d,1,1))
    noError<-F
    return(list(C=mat))
  })
  noError<-is.logical(init)
  while(!noError){
    mat<-init$C
    init<-tryCatch({
      base_chol<-chol(mat)
      noError<-T
    },error=function(e){
      print(e)
      mat<-diag(1/rgamma(d,1,1))
      return(list(C=mat))
    })
  }
  return(list(final_mat=mat,chol_mat=base_chol))
}




gen_location_general<-function(tr, sigma2, dime){
  tr$tip.label<-paste("N", 1:Ntip(tr), sep="")
  tr$node.label=paste("N",(Ntip(tr)+1):(Ntip(tr)+Nnode(tr)),sep="")
  tr_phylo4<-as(tr, "phylo4")
  location_data<-data.frame(matrix(rep(0,dime),nrow=1))
  from_vec<-rootNode(tr_phylo4)
  rownames(location_data)<-names(from_vec)
  Num_tip_located<-0
  base_chol_lt<-chol_PD(sigma2)
  sigma2<-base_chol_lt$final_mat
  base_chol<-base_chol_lt$chol_mat
  
  
  #print("finish chol decomp")
  while(length(from_vec)!=0){
    #print(from_vec)
    #print(Num_tip_located)
    next_from_nonTip<-numeric()
    for(cur_from_idx in 1:length(from_vec)){
      from_<-from_vec[cur_from_idx]
      to_<-descendants(tr_phylo4, node=names(from_), type="children")
      edge_<-edgeLength(tr_phylo4)[getEdge(tr_phylo4, node=names(to_), type="descendant")]
      from_loc<-as.matrix(location_data[names(from_),],nrow=length(from_))
      to_loc<- t(t(matrix(rnorm(dime*length(edge_),sd=sqrt(edge_)),ncol=dime) %*% base_chol) + matrix(rep(from_loc,length(edge_)),nrow=dime))
      rownames(to_loc)<-names(to_)
      colnames(to_loc)<-colnames(location_data)
      location_data<-rbind(location_data, to_loc)
      next_from_nonTip<-c(next_from_nonTip,to_[!(names(to_) %in% tipLabels(tr_phylo4))])
      Num_tip_located<-Num_tip_located+sum((names(to_) %in% tipLabels(tr_phylo4)))
    }
    from_vec<-next_from_nonTip
  }
  tr_phylo4d<-phylo4d(tr_phylo4, all.data=location_data)
  tr_phylo4d<-addData(tr_phylo4d,tip.data=t(matrix(rep(sigma2[upper.tri(sigma2,diag=T)],N),nrow=((1+d)*d/2) )))
  return(tr_phylo4d)
}

gen_location_homeskd<-function(tr, sigma2, dime){
  tr$tip.label<-paste("N", 1:Ntip(tr), sep="")
  tr$node.label=paste("N",(Ntip(tr)+1):(Ntip(tr)+Nnode(tr)),sep="")
  tr_phylo4<-as(tr, "phylo4")
  location_data<-data.frame(matrix(rep(0,dime),nrow=1))
  from_vec<-rootNode(tr_phylo4)
  rownames(location_data)<-names(from_vec)
  Num_tip_located<-0
  base_chol_lt<-chol_PD(sigma2)
  sigma2<-base_chol_lt$final_mat
  base_chol<-base_chol_lt$chol_mat
  
  
  #print("finish chol decomp")
  while(length(from_vec)!=0){
    #print(from_vec)
    #print(Num_tip_located)
    next_from_nonTip<-numeric()
    for(cur_from_idx in 1:length(from_vec)){
      from_<-from_vec[cur_from_idx]
      to_<-descendants(tr_phylo4, node=names(from_), type="children")
      edge_<-edgeLength(tr_phylo4)[getEdge(tr_phylo4, node=names(to_), type="descendant")]
      from_loc<-as.matrix(location_data[names(from_),],nrow=length(from_))
      to_loc<- t(t(matrix(rnorm(dime*length(edge_),sd=sqrt(edge_)),ncol=dime) %*% base_chol) + matrix(rep(from_loc,length(edge_)),nrow=dime))
      rownames(to_loc)<-names(to_)
      colnames(to_loc)<-colnames(location_data)
      location_data<-rbind(location_data, to_loc)
      next_from_nonTip<-c(next_from_nonTip,to_[!(names(to_) %in% tipLabels(tr_phylo4))])
      Num_tip_located<-Num_tip_located+sum((names(to_) %in% tipLabels(tr_phylo4)))
    }
    from_vec<-next_from_nonTip
  }
  tr_phylo4d<-phylo4d(tr_phylo4, all.data=location_data)
  param_df<-t(matrix(rep(sigma2[1,1],N),nrow=1 ))
  colnames(param_df)<-c("sigma2")
  tr_phylo4d<-addData(tr_phylo4d,tip.data=param_df)
  return(tr_phylo4d)
}