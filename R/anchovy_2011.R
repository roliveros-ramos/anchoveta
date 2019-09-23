

# Von Bertalanffy

VB = function(L,L.inf,k) {
  # Brody's recurrence equation for von Bertalanffy (VB) model
  # For seasonal VB use k=f(k0,WP,c)
  out=(L.inf-L)*(1-exp(-k))+L
  return(out)
}

# Natural mortality

mortality.discrete = function(L,par.m) {
  # Calculate natural mortality for discrete model give a length
  # m is a set of factors of increase, Mi=M0*prod(mi)
  # Lc is a set of cut lengths for increases in mortality
  # Example: for m=c(m1,m2,m3), Lc=c(4,8,12) ans Ma (adults)
  # M(12+)= Ma, M(8-12)=Ma*m1, M(4-8)=Ma*m1*m2 and M(4-)=Ma*m1*m2*m3
  with(par.m, {	
    m=c(1,m)
    M=Ma*cumprod(m)
    Lc=rev(c(Lc,Inf))
    output=L
    for(i in 1:length(Lc)) output[L<Lc[i]]=M[i] 
    output[L<0]=0
    return(output)
  } )
}


mortality.continuous = function(L,par.m) {
  with(par.m, {	
    output = Ma+(Mj-Ma)*(1-exp(-beta*L^(-gamma)))
    output[L<0]=0
    return(output)
  } )
}



mortality.constant = function(L,par.m) {
  with(par.m, {	
    output = rep(Ma,len=length(L))
    output[L<0]=0
    return(output)
  }
  )
}


# length-weight relationship
peso = function(L,a,b=3,unit=1E-06) {
  # unit is the convertion factor grams to tonnes
  out=unit*a*L^b
  return(out)
}

# average length of a cohort
average.length = function(N,size) {
  if(sum(N,na.rm=TRUE)==0) out=0 else out=weighted.mean(size,N)
  return(out)
}

# average sd of a cohort
average.sd = function(N,size,sd.min=0.25) {
  if(sum(N,na.rm=TRUE)==0) {
    out=0
  } else {
    ex2=weighted.mean(size^2,N)
    ex=weighted.mean(size,N)
    out=sqrt(max(ex2-ex*ex,0))
    out=max(out,sd.min)
  }
  return(out)
}

# calculate proportions
normalize = function(x) {
  suma=if(sum(x,na.rm=TRUE)!=0) sum(x) else 1
  out=x/suma
  return(out)
}

# Simple trapezoidal numerical integration
trapz=function(FUN,lower,upper,...,subdivisions=max(2,(upper-lower)/30)) {
  
  x=seq(from=lower,to=upper,by=(upper-lower)/subdivisions)
  y=FUN(x[-length(x)],...)+4*FUN(0.5*(x[-length(x)]+x[-1]),...)+FUN(x[-1],...)
  out=(1/6)*sum(diff(x)*y)
  
  return(out)
}


#### Normal density

ddnorm1=function(x,n,l,sd) {
  out=n*(1/(sqrt(2*pi)*sd))*exp(-0.5*((x-l)/sd)^2)
  return(out)
}


ddnorm2=function(x,n,l,sd) {
  out=n*dnorm(x,mean=l,sd=sd)
  return(out)
}


# cohort data to length frequencies (to generalize for any cohort dispersion) 

cohort.to.length = function(cohort,lower,upper,delta=0.5,delta.c=0.1,FUN=ddnorm1) {
  # Return a vector with length distribution of a 'cohort' (N,L,S)
  # 'lower' and 'upper' are the (class mark) limits of integration 
  # 'delta' is the length of the class interval (as measured by fishery)
  # 'delta.c' is the length for integration of the length frequency distribution
  # 'FUN' gives the density distribution (normal for one cohort)
  if(length(cohort)!=3) stop("No cohort data recognized (N,L,SD)")
  FUN	= match.fun(FUN)
  n.step	= ceiling(delta/delta.c)
  lower	= lower-0.5*delta 
  upper	= upper+0.5*delta
  n=cohort[1]; l=cohort[2]; sd=cohort[3]
  x	= seq(from=lower,to=upper,by=delta/n.step)
  y	= (1/6)*(FUN(x[-length(x)],n,l,sd)+4*FUN(0.5*(x[-length(x)]+x[-1]),n,l,sd)+FUN(x[-1],n,l,sd))
  #	y	= (1/2)*(FUN(x[-length(x)],n,l,sd)+FUN(x[-1],n,l,sd))
  out	= apply(matrix((delta/n.step)*y,nrow=n.step),2,sum)
  return(out)
}


#### population data to length frequencies
pop.to.length = function(cohort,lower,upper,delta=0.5,delta.c=0.1) {
  # Return a vector with length distribution of a population (set of 'cohort's (N,L,S))
  # 'lower' and 'upper' are the (class mark) limits of integration 
  # 'delta' is the length of the class interval (as measured by fishery)
  # 'delta.c' is the length for integration of the length frequency distribution
  if(dim(cohort)[2]!=3) stop("No cohort type data recognized (N,L,SD). Try transposing 'cohort'.")
  out=apply(cohort,1,cohort.to.length,lower,upper,delta,delta.c)
  out=apply(out,1,sum)
}


#################################################################################################
## Regime shift functions
#################################################################################################
shift.pos=function(T,reg,pos,len) {
  # 'T' is the time horizon, 'reg' is the number of regimes, 'pos' the position of the regime shift start point 
  # and 'len' is the duration of each regime transition.  
  if(length(pos)<(reg-1)) stop("Not enough regime shift start points")
  if(length(len)<(reg-1)) stop("Not enough regime shift lengths")
  out=matrix(0,ncol=reg,nrow=T)
  out[,1]=1
  for(i in 1:(reg-1)) {
    A=pmax(seq(from=1-1/(len[i]+1),by=-1/(len[i]+1),len=T-pos[i]+1),0)
    out[pos[i]:T,i]=A
    out[pos[i]:T,i+1]=1-A
  }
  return(out)
}


#################################################################################################
## Time index generation
#################################################################################################
meses=function(mes.inicial,T) {
  times=c(mes.inicial:12,rep(1:12,len=T-length(mes.inicial:12)))
  return(times)
}

anos=function(mes.inicial,year.inicial,T) {
  x=rep(year.inicial,len=length(mes.inicial:12))
  year.final=year.inicial+ceiling(T/12)
  year.inicial=year.inicial+1
  times=c(x,rep(year.inicial:year.final,each=12, len=T-length(mes.inicial:12)))
  return(times)
}

anos.cal=function(mes.inicial,years) {
  out=if(mes.inicial==1) unique(years) else unique(years)[-1] 
  return(out)
}

edad0=function(mes,rec) {
  out=mes-rec
  out[out<0]=out[out<0]+12
  out=min(out)%%12
}


#################################################################################################

