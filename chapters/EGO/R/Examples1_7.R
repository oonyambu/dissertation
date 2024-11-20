## R code for JASA paper "A Projection Space-Filling Criterion and Related Optimality Results"
## R code for Examples 1-7
## (authors) 3/11/2023

source("stratification_pattern.R") # SP() and D_SP()

regular.ha=function(k)
{# generate Hadamard matrix of order 2^k using the k-fold Kronecker product,
 # which is later used in the construction of SOAs in Examples 2-7.
  h0=matrix(c(1,1,1,-1),ncol=2,nrow=2)
  h=h0
  for( i in 1:(k-1)){
    h=kronecker(h0,h)
  }
  return(as.matrix(h))
}


#################################
#######    Example 1      #######
#################################
x1=c(0,2,3,1,6,4,5,7)
x2=c(0,6,2,4,3,5,1,7)
d=cbind(x1,x2) # GSOA(8,2,8,3) in Example 1

## making plots
par(mfrow = c(1,3),
    oma = c(1,1,0,1) + 0.1,
    mar = c(0,0.5,1,0) + 0.1)

par(pty="s",mgp=c(0,0.1,0))
x=d+0.5

# 1D plot    
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n", ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
position=seq(.5,8.5,2)
p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
polygon(c(0,8,8,0),c(0,0,1,1),col="#DCDCDC",border = NA)
polygon(c(0,8,8,0),c(2,2,3,3),col="#DCDCDC",border = NA)
polygon(c(0,8,8,0),c(4,4,5,5),col="#DCDCDC",border = NA)
polygon(c(0,8,8,0),c(6,6,7,7),col="#DCDCDC",border = NA)
abline(h =c(0,8) ,lwd = 1,lty = 3,col="gray0",lend=1)
sapply(0:8, function(i) abline(v =i ,lwd = 1,lty = 3,col="gray0",lend=1))
points(x[,1],x[,2],pch = 20,cex=0.5)

# 2D plot: 2x4  
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n", ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
sapply(0:4, function(i) abline(h =i*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
sapply(0:4, function(i) abline(v =i*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
position=seq(.5,8.5,2)
p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
points(x[,1],x[,2],pch = 20,cex=0.5)

# 2D plot: 4x2  
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n", ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
sapply(0:4, function(i) abline(h =i*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
sapply(0:4, function(i) abline(v =i*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
position=seq(.5,8.5,2)
p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
points(x[,1],x[,2],pch = 20,cex=0.5)


#################################
#######    Example 2      #######
#################################
k=3
d.full=regular.ha(k)[,-1] # generate full 2^3 factorial

a.in=c(1,2,3);b.in=c(4,4,4);c.in=c(2,1,6)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d1=2*a+b+c/2+7/2 # The first GSOA for Example 2
maxi=3*3 # the maximum i which = m*p
D_SP(d1,3,maxi) # find the stratification patterns & space-filling patterns

a.in=c(1,2,4);b.in=c(2,5,3);c.in=c(4,1,2)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d2=2*a+b+c/2+7/2 # The second GSOA for Example 2
D_SP(d2,3,maxi) # find the stratification patterns & space-filling patterns

## making plots
D=d1 # or D=d2
x=D+0.5

dev.new(width=10, height=4)
par(mfrow = c(2,3),
    oma = c(1,1,0,1) + 0.1,
    mar = c(0,0.5,1,0) + 0.1)

par(pty="s",mgp=c(0,0.1,0))

# 2D:2x4    
for (i in 1:2){
  for(j in (i+1):3){
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n", ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
sapply(0:4, function(i) abline(h =i*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
sapply(0:4, function(i) abline(v = i*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
position=seq(.5,8.5,2); p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
points(x[,i],x[,j],pch = 20,cex=0.5)
}
}

# 2D:4x2 
for (i in 1:2){
  for(j in (i+1):3){
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n", ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
sapply(0:4, function(i) abline(h =i*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
sapply(0:4, function(i) abline(v = i*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
position=seq(.5,8.5,2);p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
points(x[,i],x[,j],pch = 20,cex=0.5)
  }
}

# the shadow in 4x2 plot of pair(1,2) for D2
polygon(c(0,2,2,0),c(4,4,8,8),col="#DCDCDC",lty=3,lend=1)
polygon(c(2,4,4,2),c(0,0,4,4),col="#DCDCDC",lty=3,lend=1)
polygon(c(4,6,6,4),c(4,4,8,8),col="#DCDCDC",lty=3,lend=1)
polygon(c(6,8,8,6),c(0,0,4,4),col="#DCDCDC",lty=3,lend=1)


#################################
#######    Example 3      #######
#################################
k=5
d.full=regular.ha(k)[,-1] # generate full 2^k factorial

a.in=c(27,25,24,28);b.in=c(9,15,9,12);c.in=c(31,31,31,24)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d3=2*a+b+c/2+7/2 
sapply(2:4,function(j) SP(d3,3,4,j)) # find (P_42,P_43,P_44) for d3
D_SP(d3,3,6)[[2]] # find S1-S6 for d3

a.in=c(1,2,4,8,16,15,19,21,22);b.in=c(24,28,3,5,10,6,12,14,11);c.in=c(6,5,10,20,7,12,7,18,17)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d=2*a+b+c/2+7/2 ## SOA(32,9,8,3) in Table 6
d4=d[,c(5,6,7,8)];d4
sapply(2:4,function(j) SP(d4,3,4,j)) # find (P_42,P_43,P_44) for d4
D_SP(d4,3,6)[[2]] # find S1-S6 for d4

## make plots
d=d3 # or d=d4

dev.new(width=10, height=4)
par(mfrow = c(3,6),
    oma = c(1,1,0,1) + 0.1,
    mar = c(0,0.5,1,0) + 0.1)

par(pty="s",mgp=c(0,0.1,0))
x=d+0.5

# 2D plot:2x8
for (i in 1:3){
  for(j in (i+1):4){
plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n",ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
sapply(0:8, function(l) abline(h =l*1 ,lwd = 1,lty = 3,col="gray0",lend=1))
sapply(0:2, function(l) abline(v =l*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
position=seq(.5,8.5,2)
p.label=seq(0,8,2)
axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
points(x[,i],x[,j],pch = 20,cex=0.5)
  }
}

# 2D plot:8x2
for (i in 1:3){
  for(j in (i+1):4){
    plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n",ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
    sapply(0:2, function(l) abline(h =l*4 ,lwd = 1,lty = 3,col="gray0",lend=1))
    sapply(0:8, function(l) abline(v =l*1 ,lwd = 1,lty = 3,col="gray0",lend=1))
    position=seq(.5,8.5,2)
    p.label=seq(0,8,2)
    axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
    axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
    points(x[,i],x[,j],pch = 20,cex=0.5)
  }
}

# 2D plot:4x4
for (i in 1:3){
  for(j in (i+1):4){
    plot(0:8,lwd=2,col="black",type="n",xaxt="n",yaxt="n",ylim=c(-0.5,8.5),xlim=c(-0.5,8.5),ann=FALSE)
    sapply(0:4, function(l) abline(h =l*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
    sapply(0:4, function(l) abline(v =l*2 ,lwd = 1,lty = 3,col="gray0",lend=1))
    position=seq(.5,8.5,2)
    p.label=seq(0,8,2)
    axis(1,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
    axis(2,at=position,labels=p.label,tck=-.03,cex.axis=0.6,lwd.ticks=0.5)
    points(x[,i],x[,j],pch = 20,cex=0.5)
  }
}

# the shadow in 2x8 plot of pair(3,4) for D3
polygon(c(0,4,4,0),c(1,1,2,2),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(0,4,4,0),c(3,3,4,4),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(0,4,4,0),c(5,5,6,6),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(0,4,4,0),c(7,7,8,8),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,8,8,4),c(0,0,1,1),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,8,8,4),c(2,2,3,3),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,8,8,4),c(4,4,5,5),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,8,8,4),c(6,6,7,7),col="#DCDCDC",border =  "gray0",lty=3,lend=1)

# the shadow in 4x4 plot of pair(1,3) for D3
polygon(c(0,2,2,0),c(2,2,4,4),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(0,2,2,0),c(6,6,8,8),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(2,4,4,2),c(0,0,2,2),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(2,4,4,2),c(4,4,6,6),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,6,6,4),c(2,2,4,4),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(4,6,6,4),c(6,6,8,8),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(6,8,8,6),c(0,0,2,2),col="#DCDCDC",border =  "gray0",lty=3,lend=1)
polygon(c(6,8,8,6),c(4,4,6,6),col="#DCDCDC",border =  "gray0",lty=3,lend=1)

#################################
#######    Example 4      #######
#################################
k=5
d.full=regular.ha(k)[,-1]
a.in=c(1,2,4,8,16,7,19,29);b.in=c(24,20,9,6,5,27,12,3);c.in=c(11,8,14,2,26,13,25,23)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d=2*a+b+c/2+7/2 #SOA(32,2^3,8,3)

list.subd=function(d,p,m)
{# rank all subdesigns of m factors according to their stratification patterns & space-filling patterns, respectively
  M=ncol(d)
  imax=p*m # the maximum value for i
  com=combn(M,m)
  d.sub=d[,com[,1]]
  result=D_SP(d.sub,p,imax)
  r1=result[[1]]
  r2=result[[2]]
  for(l in 2:ncol(com)){
    res=D_SP(d[,com[,l]],p,imax)
    r1=rbind(r1,res[[1]][3,])
    r2=rbind(r2,res[[2]][2,])
  }
  rr1=cbind(r1[-c(1,2),],NA,t(com))
  rrr1=rr1[do.call(order, c(decreasing = FALSE, data.frame(rr1[,1:ncol(r1)]))),]
  rrr1=rbind(cbind(r1[c(1,2),],matrix(NA,2,(m+1))),rrr1)
  rr2=cbind(r2[-1,],NA,t(com))
  rrr2=rr2[do.call(order, c(decreasing = FALSE, data.frame(rr2[,1:ncol(r2)]))),]
  rrr2=rbind(c(r2[1,],rep(NA,(m+1))),rrr2)
  return(list(rrr1,rrr2)) # the last m+1 columns of rrr1 & rrr2 are the column indices in subdesigns  
  }

m=3
subdesigns.list=sapply(1:m,function(l) list.subd(d,3,l))
D_SP(d,3,24) #m=8

#################################
#######    Example 5      #######
#################################
k=5
d.full=regular.ha(k)[,-1]

a.in=c(1,2,4,8,16,7,11,19,29);b.in=c(24,20,9,6,5,27,17,12,3);c.in=c(rep(10,9))
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d5=2*a+b+c/2+7/2 #SOA(32,9,8,3) in Example 1 of Shi and Tang(2020)
sapply(2:4,function(j) SP(d5,3,4,j)) # find (P_42,P_43,P_44) 


a.in=c(1,2,4,8,16,15,19,21,22);b.in=c(24,28,3,5,10,6,12,14,11);c.in=c(6,5,10,20,7,12,7,18,17)
a=d.full[,a.in];b=d.full[,b.in];c=d.full[,c.in]
d6=2*a+b+c/2+7/2 #SOA(32,9,8,3) in Table 6
sapply(2:4,function(j) SP(d6,3,4,j)) # find (P_42,P_43,P_44) 

#################################
#######    Example 6      #######
#################################
m=4
## d5:SOA(32,9,8,3) in Example 1 of Shi and Tang(2020)
subdesigns.list=sapply(1:m,function(l) list.subd(d5,3,l)) # list.subd() is given under Example 4
D_SP(d5,3,27) 

## d6:SOA(32,9,8,3) in Table 6
subdesigns.list=sapply(1:m,function(l) list.subd(d6,3,l)) # list.subd() is given under Example 4
D_SP(d6,3,27) 

#################################
#######    Example 7      #######
#################################
k=4
d.full=regular.ha(k)[,-1]
#############################
a.in=9:15;b.in=rep(8,7)
c1.in=c(2:7,1);c2.in=c(2,4,2,2,1,1,1);c3.in=c(10,rep(9,6))  
a=d.full[,a.in];b=d.full[,b.in]
c1=d.full[,c1.in];c2=d.full[,c2.in];c3=d.full[,c3.in]
d1=2*a+b+c1/2+7/2;d2=2*a+b+c2/2+7/2;d3=2*a+b+c3/2+7/2
D_SP(d1,3,6)[[1]];D_SP(d2,3,6)[[1]];D_SP(d3,3,6)[[1]]

