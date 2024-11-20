## R code for JASA paper "A Projection Space-Filling Criterion and Related Optimality Results"
## main function: SP() and D_SP() called by Examples1_7.R 
## (authors) 3/11/2023

## functions for calculating P_ij's and S_i's for GSOAs with s=2 and p=3
x1=rep(c(1,-1),each = 4);
x2=rep(c(1,-1),2,each = 2);
x3=rep(c(1,-1),4);
chi=as.matrix(cbind(rep(1,8),x1,x2,x1*x2,x3, x1*x3, x2*x3, x1*x2*x3))
colnames(chi) <- NULL #chi includes characters for any x, u in {0,...,2^p} with p=3

weights=function(i,j,p)
{ # find all possible weights for \rho=i(subregion size) & \psi=j(dimension)
  res=gtools::permutations(n=p,r=j,repeats.allowed=T)
  wg=matrix(res[which(rowSums(res)==i),],nrow=length(which(rowSums(res)==i)))
  return(wg)
}

mapping.level=function(weight)
{ # map the weights to levels of the design
  x=sapply(weight,function(l) 2^(l-1):(2^l-1)) # 2 is s.
  if(is.list(x)){ value=expand.grid(x)
  }else if(is.matrix(x)){value=expand.grid(data.frame(x))
  }else{value=weight }
  colnames(value)<-NULL
  return(data.frame(value))
}

sp=function(d,chi,u) # u is the weights
{ # This function finds the stratification pattern for d and j=ncol(d) and is used in SP().
  # u is the weights corresponding i&j, where i is given by the arguments of SP().
  res=rep(0,nrow(u))
  for(l in 1:nrow(u)){
    su=1
    for(j in 1:ncol(u)){
      su=su*chi[(d[,j]+1),(u[l,j]+1)]
    }
    res[l]=abs(sum(su)/nrow(d))
  }
  return(res)
}

SP=function(D,p,i,j)
{ # find the stratification pattern P_ij for any D
  D=as.matrix(D)
  m=ncol(D)
  if(j>i && j<=m){ # not available values
    return("NA")
  }else if(i>p && j<ceiling(i/p)){ # not available values
    return("NA")
  }else if(j>m){ # not defined values
    return("-")
  }else{
    wei=weights(i,j,p)
    y=apply(wei,1,mapping.level)
    wt=rlist::list.rbind(y)
    if(dim(wt)[2]==1 && length(unique(wt[,1]))==1){wt=t(wt)}
    ind=combn(m,j)
    s=apply(ind,2, function(l) sp(as.matrix(D[,l]),chi,wt))
    return(sum(s))
  }
}

D_SP=function(D,p,maxi)
{ # find P_ij and S_i for i=1,...,maxi and j=1,...,i.
  # the P_ij's that are missing are omitted.
  res=matrix(c(1,1,0),nrow=3) 
  res2=matrix(c(1,0),nrow=2)  
  for(i in 2:maxi){
    ind=matrix(c(rep(i,i),1:i),nrow=2,byrow=TRUE)
    x=sapply(1:i,function(l) SP(D,p,i,l))
    ind.new=as.matrix(ind[,complete.cases(suppressWarnings(as.numeric(x)))])
    res=cbind(res,rbind(ind.new,na.omit(suppressWarnings(as.numeric(x)))))
    res2=cbind(res2,c(i,sum(na.omit(suppressWarnings(as.numeric(x))))))
  }
  rownames(res) <- c("i", "j","P_ij")
  rownames(res2) <- c("i", "S_i")
  return(list(res,res2))
}
 

