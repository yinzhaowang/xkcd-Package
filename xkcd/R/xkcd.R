

################### rxkcd 
rxkcd=function(n, sd=1){
    mean = 0
    x=rnorm(n,mean,sd)
    temp=dnorm(x,mean,sd)
    y=runif(n,mean,temp)
    return(y)
}

################### dxcd
dxkcd=function(x, sd = 1, log = FALSE, swap.end.points = FALSE){
    # some checks
    stopifnot(sum(is.na(x))+sum(is.na(sd))==0)
    stopifnot(sum(sd<=0)==0 )
    stopifnot(sum(is.infinite(sd))==0)
        
    if(length(sd) == 1){
        sd = rep(sd, length(x))
    }
    if(length(sd) < length(x)){
        sd = append(sd, rep(1, length(length(x) - length(sd))))
    }
    sq2psd=sqrt(2*pi)*sd
    up_bd=1/sq2psd
    log2sd2=log(2*sd^2) # const2
    log2psd2=log2sd2+log(pi)
    y=numeric(length(x))
    

    
     
    ########## consider the case when input x is inside domain
    i= (x>0) & (x<up_bd)
    if (sum(i)>0){
        if (log == FALSE & swap.end.points == FALSE){
            
            y[i]=2*sqrt((-1/2*log2psd2[i]-log(x[i]))*2*sd[i]^2)
        }
        else if(log == FALSE & swap.end.points == TRUE){
            
            y[i]=2*sqrt( -2*sd^2*log1p(-sq2psd[i]*x[i])  )
        }
        else if(log == TRUE & swap.end.points == FALSE){
            
            y[i]=log(2)+1/2*( log(-1/2*log2sd2[i]-1/2*log(pi)-log(x[i]))  +log2sd2[i])
        }
        else if(log == TRUE & swap.end.points == TRUE){
            
            y[i]=log(2)+1/2*(log( -log1p(-sq2psd[i]*x[i]))+log2sd2[i] )
        }
    }  
    
    #avoid extra calculation
    i_0= (x==0)
    if (sum(i_0)>0 ){
        if (log == FALSE & swap.end.points == FALSE){y[i_0]=Inf}
        else if (log == FALSE & swap.end.points == TRUE){y[i_0]=0}
        else if (log == TRUE & swap.end.points == FALSE){y[i_0]=Inf}
        else if (log == TRUE & swap.end.points == TRUE){y[i_0]=-Inf}
    }
    i_up= (x==up_bd)
    if (sum(i_up)>0 ){
        if (log == FALSE & swap.end.points == FALSE){y[i_up]=0}
        else if (log == FALSE & swap.end.points == TRUE){y[i_up]=Inf}
        else if (log == TRUE & swap.end.points == FALSE){y[i_up]=-Inf}
        else if (log == TRUE & swap.end.points == TRUE){y[i_up]=-Inf}
    }
    
    return(y)
}

################ pxkcd 

pxkcd=function(q, sd = 1, log.p = FALSE, swap.end.points = FALSE){
    stopifnot(sum(is.na(1))+sum(is.na(sd))==0)
    stopifnot(sum(sd<=0)==0 )
    stopifnot(sum(is.infinite(sd))==0)
    if(length(sd) == 1){
        sd = rep(sd, length(q))
    }
    if(length(sd) < length(q)){
        sd = append(sd, rep(1, length(length(q) - length(sd))))
    }
    

    p=numeric(length(q))
    sq2psd=sqrt(2*pi)*sd
    up_bd=1/sq2psd
    log2sd2=log(2*sd^2)
    log2psd2=log2sd2+log(pi)

    
    # check boundaries 
    
    in0= (q>0) & (q<up_bd)
    i2= (q>0)&(q<2^{-30})
    i3=(q>=2^{-30})&(q< up_bd)
    
    i_over = (q>=up_bd)
    p[i_over] = 1
    
    if (sum(in0)>0){

          
        #case1
        if (log.p == FALSE & swap.end.points == FALSE){
            temp1=dxkcd(q[in0],sd=sd[in0])/2 
            temp2=pnorm(-temp1 ,mean=0,sd=sd[in0]) 
            p[in0]=2*temp2+q[in0]*(2)*temp1
        }
        #case2
        else if(log.p == FALSE & swap.end.points == TRUE){
            if(sum(i2)>0){
                c=sqrt(2*sd^3 *sqrt(2*pi) )*4/3
                p[i2]=c*q[i2]^{3/2}
            }
            if(sum(i3)>0){
                temp11=dxkcd(q[i3],sd=sd[i3],swap.end.points = TRUE)/2 
                temp22=pnorm(-temp11 ,mean=0,sd=sd[i3])
                p[i3]=1-2*(temp22+(up_bd[i3]-q[i3])*temp11)
            }
        }
        #case3
        else if(log.p == TRUE & swap.end.points == FALSE){
            i4= (q>0)&(q<2^{-1000})
            i5=(q>=2^{-1000})&(q< up_bd)

            if(sum(i4)>0){
                temp3=log(sqrt(2* pi)*sd[i4] *q[i4] )
                p[i4]=log(2/sqrt(pi))+ temp3 + log(1/2-temp3)-1/2*log(-temp3)
            }
            if(sum(i5)>0){
                temp111=dxkcd(q[i5],sd=sd[i5])/2
                temp222=pnorm(-temp111 ,mean=0,sd=sd[i5])
                t=2*(temp222+q[i5]*temp111)
                p[i5]=log(t)
            }      
        }
        #case4
        else if(log.p == TRUE & swap.end.points == TRUE){
            if(sum(i2)>0){
                p[i2]=log(sqrt(2*sd[i2]^3 *sqrt(2*pi) )*4/3)+3/2*log(q[i2])  
            }
            if(sum(i3)>0){
                temp11=dxkcd(q[i3],sd=sd[i3],swap.end.points = TRUE)/2 
                temp22=pnorm(-temp11 ,mean=0,sd=sd)
                p[i3]=log(1-2*(temp22+(up_bd[i3]-q[i3])*temp11))
            }
        }
    }   
    return(p)
}


##############  qxkcd


qxkcd=function(p, sd = 1, log.p = FALSE,swap.end.points = FALSE){
    sq2psd=sqrt(2*pi)*sd
    up_bd=1/sq2psd
    
    stopifnot(sum(is.na(p))+sum(is.na(sd))==0)
    stopifnot(sum(sd<=0)==0 )
    stopifnot(sum(is.infinite(sd))==0)
    stopifnot(sum(p<0)==0 | log.p==FALSE)
    
    q=numeric(length(p))
    for(i in 1:length(p)){
            fun=function(q){return(pxkcd(q, sd = sd, log.p = log.p,swap.end.points = swap.end.points)- p[i])}
            q[i]=uniroot(fun, lower=0, upper=up_bd)$root
        }
    return(q)
    
}


