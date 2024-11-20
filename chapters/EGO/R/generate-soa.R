## R code for JASA paper "A Projection Space-Filling Criterion and Related Optimality Results"
## main function: get.soa.new() called by km-simuation.R 
## (authors) 3/11/2023

## functions to generate SOAs
ffd2gen = function(gen)
{ # generate 2-level ffds with generator matrix 
    if(is.vector(gen) || nrow(gen) == 1)     rbind(0, gen) %% 2
    else{
        g <- gen[nrow(gen),] # last row
        x <- ffd2gen(gen[-nrow(gen),])  
        x <- t(x)
        x<- cbind(x, x+g ) %% 2
        x<- t(x)
        dimnames(x) <- list(1:nrow(x), 1:ncol(x))
        x
        }
}

ffd.s = function(n, s=2, rev=T)
{ # a simple function to generate 2-level ffd
	if(s !=2) stop("s !=2: not implemented")
	Identity = function(n) diag(rep(1,n))
    gen <- t(ffd2gen(Identity(n)))[,-1]	# remove col of 0's
    if(rev) x <- ffd2gen(gen[n:1,])  # reverse the order so that it is easy to read
    else x <- ffd2gen(gen)  
	list(code=x)
}

get.xyz=function(k)
{ # create columns of X, Y, Z using Yates orders
	if(k==2) return(cbind(x=c(1,2,3), y=c(2,3,1), z=c(3,1,2)))
	else if(k==3) return(cbind(x=c(1:7), y=c(7,5,2,1,6,4,3), z=c(6,7,1,5,3,2,4)))
	xyz=get.xyz(k-2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	x=c(0, xyz[,1]); x=c(x, k1+x, k2+x, k3+x)[-1]
	y=c(0, xyz[,2]); y=c(y, k2+y, k3+y, k1+y)[-1]
	z=c(0, xyz[,3]); z=c(z, k3+z, k1+z, k2+z)[-1]
	cbind(x=x,y=y,z=z)
}
soa.ST5=function(k, iC=k1, printC=F)
{ # construct SOA(n=2^k,n/4-1,8, 3) as in Shi and Tang (2020, Theorem 5)
	# iC=k1=n/4 is the unique choice.
	xyz = get.xyz(k-2)
	a=ffd.s(k,2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	A=a$code[,c(k1+xyz[,1])]
	B=a$code[,c(k2+xyz[,2])]
	B1=a$code[,c(k3+xyz[,3] )] # B1 =A+B
	setAB=c(k1+xyz[,1], k2+xyz[,2], k3+xyz[,3], xyz[,1]) # A, B, B1, A2
	setC=(1:(2^k-1))[-setAB] # all possible C not in (A, A2, B, B1)
	if(printC) print(setC)
#	cc=rep(setC, ncol(A))[1:ncol(A)]
	if(missing(iC)) iC=setC[1]
	cc=rep(iC, ncol(A))[1:ncol(A)]
	C=a$code[, cc ] # avoid C in 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

soa.ST5x=function(k, iC=1, printC=F)
{ # construct SOA(n=2^k,n/4,8, 3) as in Shi and Tang (2020, corollary 1)
	# iC can be from 1:(k1-1)
	xyz = get.xyz(k-2)
	a=ffd.s(k,2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	A=a$code[,c(k1, k1+xyz[,1])]
	B=a$code[,c(k2, k2+xyz[,2])]
	B1=a$code[,c(k3, k3+xyz[,3] )] # B1 =A+B
	setAB=c(k1+xyz[,1], k2+xyz[,2], k3+xyz[,3], xyz[,1]) # A, B, B1, A2
#	setC=(1:(2^k-1))[-setAB] # all possible C not in (A, A2, B, B1)
#	if(printC) print(setC)
#	cc=rep(setC, ncol(A))[1:ncol(A)]
#	if(missing(iC)) iC=setC[1]
	cc=rep(iC, ncol(A))[1:ncol(A)]
	C=a$code[, cc ] # avoid C in 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

get.abc=function(k)
{ # initial A, B, C as in Section 4.1 for constructing SOA(n=2^k,5n/16,8, 3)	 or SOA(32,9,8,3)		
	if(k==4) return(cbind(x=c(1,2,4,8,15), y=c(12,9,3,6,5), z=c(11,7,9,5,12))) 
	else if(k==5) return(cbind(x=c(1,2,4,8,16,15,19,21,22), y=c(24,28,3,5,10,6,12,14,11), z=c(6,5,10,20,7,12,7,18,17)))
	else if(k==7){
	# create columns of A, B, C using Yates orders as Theorem 3 for N=128
A=c(1,2,4,8,15,17,18,20,24,31,33,34,36,40,47,49,50,52,56,63,65,66,68,72,79,81,82,84,88,95,97,98,100,104,111,113,114,116,120,127);
					B=c(42,37,25,3,117,74,41,10,14,102,92,69,23,6,83,90,73,71,21,86,54,28,7,5,57,61,44,26,19,53,60,12,9,13,58,55,62,35,27,38);
					C=c(60,27,90,69,54,61,119,71,59,108,75,105,77,85,53,76,78,125,59,106,61,101,122,107,53,42,119,109,46,102,74,9,3,6,73,10,69,42,5,21);
	 return(cbind(x=A, y=B, z=C))
	}

	if(k<4) stop("k<4")
	xyz=get.abc(k-2)
	k0=0; k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	x=xyz[,1]; x=c(x, k1+x, k2+x, k3+x)
	y=xyz[,2]; y=c(y, k2+y, k3+y, k1+y) # original as in Shi and Tang (2020)
	z=xyz[,3]; z=c(z, k2+z, k3+z, k1+z)
	cbind(x=x,y=y,z=z)
}

soa.ST3=function(k)
{ # constructing SOA(n=2^k,5n/16,8, 3) as Theorem 3 in Section 4.1
# except for k=5, soa(32,9,8,3)
	xyz = get.abc(k)
	a=ffd.s(k,2)
	A=a$code[,(xyz[,1])]
	B=a$code[,(xyz[,2])]
	C=a$code[,(xyz[,3])] # 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

get.soa.new=function(N=64, k=15, rand=F)
{ 
	p=log(N,2);
	if(N != 2^p) stop("N !=2^p")
	if(k>N/4) x=soa.ST3(p) # Theorem 3, m=5N/16
	else if(k==N/4) x=soa.ST5x(p) # Shi and Tang (2020, corollary 1)
	else x=soa.ST5(p) # Theorem 6, m=N/4-1; Shi and Tang (2020, Theorem 5)
	if(rand) x=x[,sample(ncol(x), k)]
	x[,1:k]	
}

## generate SOAs used in the simulations
# D15=get.soa.new(N=64,k=15) # soa64x15
# D20=get.soa.new(N=64,k=20) # soa64x20
# D31=get.soa.new(N=128,k=31) # soa128x31
# D40=get.soa.new(N=128,k=40) # soa128x40

## get P_ij for i=1,2,3,4 for these SOAs 
#source("stratification_pattern.R") # D_SP()
# D_SP(D15-1,3,4)
# D_SP(D20[,c(1:15)]-1,3,4) # use first 15 columns
# D_SP(D31-1,3,4)
# D_SP(D40[,c(1:31)]-1,3,4) # use first 31 columns


## functions used to search best space-filling designs
## 
search.maximin=function(N, k, t=1, nRep=100, write=T)
{
	require(SLHD)
  ## search maximinLHD nRep times
  print(date())
  a0=maximinSLHD(t=t,N,k); a0$measure; min(dist(a0$Design))
  print(date())
  if(nRep>0) nRep = nRep -1 
  for(i in 1:nRep){ a=maximinSLHD(1,N,k); if(a$measure < a0$measure) a0=a }
  a0$measure; min(dist(a0$Design))
  print(date())

  LD4=a0$Design
  dimnames(LD4)[[2]]=LETTERS[1:k]
  file4=paste("maximin", N, "x", k, ".txt", sep="")
  if(write) write.table(LD4,file4, qu=F, row=F)
  LD4
}

search.maxpro=function(N, k, nRep=100, write=T)
{
	require(maxpro)
  ## search  MaxProLHD nRep times  
  print(date())
  a0=MaxProLHD(n=N, p=k); a0$measure; min(dist(a0$Design))
  print(date())
  
  if(nRep>0) nRep = nRep -1 
  for(i in 1:nRep) {a=MaxProLHD(n=N, p=k); if(a$measure < a0$measure) a0=a }
  print(date())
  a0$measure; min(dist(a0$Design))
  LD5 = a0$Design*N+0.5 # make levels to 1-N
  dimnames(LD5)[[2]]=LETTERS[1:k]
  file5=paste("maxpro", N, "x", k, ".txt", sep="")
  if(write)write.table(LD5,file5, qu=F, row=F)
  LD5
}

search.UD=function(N, k, q=N, nRep=10, crit="CD2", write=T)
{
	require(UniDOE)
  ## search maximinLHD nRep times
  print(date())
  a0=GenUD(N, k, q, crit=crit); a0$criterion_value; 
  print(date())
  if(nRep>0) nRep = nRep -1 
  for(i in 1:nRep){ a=GenUD(N, k, q, crit=crit); if(a$criterion_value < a0$criterion_value) a0=a }
  a0$criterion_value; 
  print(date())

  LD4=a0$final_design
  dimnames(LD4)[[2]]=LETTERS[1:k]
  file4=paste("ud", N, "x", k, ".txt", sep="")
  if(q !=N) file4=paste("ud", q,"-", N, "x", k, ".txt", sep="")
  if(write) write.table(LD4,file4, qu=F, row=F)
  LD4
}

