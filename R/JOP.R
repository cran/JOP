.packageName<-'JOP'

  ########################################################
  ####### Transformation to polar coordinates ############
  ########################################################
  
  
  trafopar<-function(x)
  {
      n   <- length(x)
      if (n==1)
      {
      return(cat("Error! Function needs more dimensions!\n"))
      }
      
      
      y   <- numeric(n)
      r   <- x[1]
      phi <- numeric(length(n-1))
      
      for (i in 1:(n-1))
      {
          phi[i] <- x[i+1]
      }
      
      y[n] <- r * sin(phi[n-1])
      y[1] <- r * cos(phi[1])
      y[2] <- r * sin(phi[1])
      
        if (n==2)
      {
      return(y)
      }
      
      for (i in 3:n)
      {
          for (j in 1:(i-1))
          {
          y[j] <- y[j]*cos(phi[i-1])      
          }
      
      y[i] <- r * sin(phi[i-1])  
      }
      
      
      return(y)
  }
        
       ###############################################
       ###############################################
       #########       Function 'JOP'       ##########
       ###############################################
       ###############################################
    
  JOP<-function(nx=2,ny=1,Wstart=1,Wend=1,numbW=1,d=c(1,0),optreg=0,data=NULL,tau=NULL,mean.model=NULL,var.model=NULL,solver=0)
  {
        ## solver=0: nlminb
        ## solver=1: rgenoud
        ## optreg=0: Optimization region is a sphere
        ## optreg=1: Optimization region is a cube
        
        ## mean.model, var.model : List of functions of type mean.model<-function(x)
        ##                           that return a single response value
        
       if(is.null(data))
       {
        return("data required!!")
       }     
       if(is.null(mean.model[[1]]))
       {
        return("Type in the model function for Mean!")
       }
       if(is.null(var.model[[1]]))
       {
        return("Type in the model function for Dispersion!")
       }
  
       if(!is.list(mean.model) && !is.list(var.model))
       {
        return("mean.model and var.model have to be of type 'list'!")
       }
       if(length(d)!=ny)
       {
        return("Check Dimension of d! Have to agree with number of responses")
       }
       if(is.null(tau))
       {
        return("Type in the target values, please!")
       }
       if(solver!=0 && solver!=1 && solver!=2)
       {
         cat("Choose a solver! '0' for 'nlminb' or '1' for goslnp or '2' for genoud....")
       }
       ##################################
       ##################################
       #### Weight Matrices and     #####
       #### Standardization Matrix  #####
       ##################################
       ##################################
    
       #### Weight Matrices W   ####
       ####(Diagonal Matrices)  ####
       
       if (Wend < Wstart)
       { 
        return("Wstart must be smaller than Wend!")  
       }
       W<-vector("list",numbW)
       logat<-seq(Wstart,Wend,(Wend-Wstart)/(numbW-1))
       for(i in 1:(numbW))
       {
         W[[i]]<-diag(exp(d*logat[i])) 
       }
       cat("Weightmatrices build...\n")
       cat("\n")
       #### Standardization Matrix A ####
       #### (Diagonal Matrix)        ####

       #Designmatrix
       dimdaten<-dim(data)[1]
       xdesign<-matrix(NaN,ncol=nx,nrow=dimdaten)
       xdesign<-data[,1:nx]  
       diagA<-NULL
       
       zw<-NULL
       for(i in 1:ny)
       {
         for(j in 1:dimdaten)
         {
             zw[j]<-as.numeric(var.model[[i]](xdesign[j,]))  
         }
         diagA[i]<-1/sqrt(mean(zw))
         zw<-NULL 
       }  
       A<-diag(diagA)

       #Cost Matrix C
       for(i in 1:length(W))
       {
         W[[i]]<-t(A)%*%W[[i]]%*%A
       }   
       cat("...Standadization Matrix build...\n")
       cat("\n")
  
       
       #################################
       #################################
       ######### Optimization ##########
       #################################
       #################################
       
    ## Matrix with optimal Parameter Settings
    optmatrix<-matrix(NaN,ncol=nx,nrow=numbW)
  
    
    ## Matrix with optimal Responses
    reoptmatrix<-matrix(NaN,ncol=ny,nrow=numbW)
     
    ## Restrictions
    if(optreg==0)
    {
      lower<-c(0,rep(-pi,nx-1))
      upper<-c(max(abs(xdesign)),rep(pi,nx-1))
      Domain<-cbind(lower,upper)
      cat("...Optimization starts...\n")
      cat("\n")
      deviation<-matrix(NaN,ncol=ny,nrow=numbW)
      optval<-NULL
  
      for(i in 1:numbW)
      {
        riscfun<-function(x)
        {
          varval<-NULL
          meanval<-NULL
          for(j in 1:ny)
          {
             varval[j]<-as.numeric(var.model[[j]](trafopar(x)))
             meanval[j]<-as.numeric(mean.model[[j]](trafopar(x)))    
          }
          return(sum(diag(W[[i]]%*%diag(varval)))+sum(t(meanval-tau)%*%W[[i]]%*%(meanval-tau)))  
        }
        #Optimization is conducted on a sphere
    
        # User can choose a solver
        if(solver==0)
        { 
          opt1<-nlminb(c(max(xdesign)/2,rep(-pi/2,nx-1)),riscfun,lower=lower,upper=upper)
          opt<-opt1
          opt2<-nlminb(c(max(xdesign)/2,rep(0,nx-1)),riscfun,lower=lower,upper=upper)
          if(opt$objective>opt1$objective)
          {
            opt<-opt2
          }
          opt3<-nlminb(c(max(xdesign)/2,rep(pi/2,nx-1)),riscfun,lower=lower,upper=upper)    
          if(opt$objective>opt3$objective)
          {
            opt<-opt3
          }
          optval[i]<-opt$objective
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-trafopar(opt$par)
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](trafopar(opt$par)))
            deviation[i,k]<-as.numeric(var.model[[k]](trafopar(opt$par)))
          }
        }
        if(solver==1)
        { 
          opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2)
          optval[i]<-opt$values[length(opt$values)]
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-trafopar(opt$pars)
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](trafopar(opt$pars)))
            deviation[i,k]<-as.numeric(var.model[[k]](trafopar(opt$pars)))
          }
        }
        if(solver==2)
        {
          opt<-genoud(riscfun,nvars=nx,Domains=Domain, print.level=0,
                   boundary.enforcement=2,wait.generations=50)
          optval[i]<-opt$value
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-trafopar(opt$par)
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](trafopar(opt$par)))
            deviation[i,k]<-as.numeric(var.model[[k]](trafopar(opt$par)))
          }
        }
                                                     
      }
    }
    if(optreg==1)
    {
      lower<-rep(-max(abs(xdesign)),nx)
      upper<-rep(max(abs(xdesign)),nx)
      Domain<-cbind(lower,upper)
      cat("...Optimization starts...\n")
      cat("\n")
      deviation<-matrix(NaN,ncol=ny,nrow=numbW)
      optval<-NULL
  
      for(i in 1:numbW)
      {
        riscfun<-function(x)
        {
          varval<-NULL
          meanval<-NULL
          for(j in 1:ny)
          {
             varval[j]<-as.numeric(var.model[[j]](x))
             meanval[j]<-as.numeric(mean.model[[j]](x))    
          }
          return(sum(diag(W[[i]]%*%diag(varval)))+sum(t(meanval-tau)%*%W[[i]]%*%(meanval-tau)))  
        }
        #Optimization is conducted on a sphere
    
        # User can choose a solver
        if(solver==0)
        { 
          opt1<-nlminb(rep(max(xdesign)/2,nx),riscfun,lower=lower,upper=upper)
          opt<-opt1
          opt2<-nlminb(rep(-max(xdesign)/2,nx),riscfun,lower=lower,upper=upper)
          if(opt$objective>opt1$objective)
          {
            opt<-opt2
          }
          opt3<-nlminb(rep(0,nx),riscfun,lower=lower,upper=upper)    
          if(opt$objective>opt3$objective)
          {
            opt<-opt3
          }
          optval[i]<-opt$objective
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-opt$par
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](opt$par))
            deviation[i,k]<-as.numeric(var.model[[k]](opt$par))
          }
        }
        if(solver==1)
        { 
          opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2)
          optval[i]<-opt$values[length(opt$values)]
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-opt$pars
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](opt$pars))
            deviation[i,k]<-as.numeric(var.model[[k]](opt$pars))
          }
        }
        if(solver==2)
        {
          opt<-genoud(riscfun,nvars=nx,Domains=Domain, print.level=0,
                   boundary.enforcement=2,wait.generations=50)
          optval[i]<-opt$value
          # The results are stored in optmatrix, reoptmatrix and deviation
          optmatrix[i,]<-opt$par
          for(k in 1:ny)
          {
            reoptmatrix[i,k]<-as.numeric(mean.model[[k]](opt$par))
            deviation[i,k]<-as.numeric(var.model[[k]](opt$par))
          }
        }
                                                     
      }
    }
    deviation<-sqrt(deviation)
    cat("...Optimal Parameters and Responses received...\n")
    cat("\n")
 
    #################################
    #################################
    #### Joint optimization Plot ####
    #################################
    #################################
    
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(signif(min(xdesign),digits=3),signif(max(xdesign),digits=3),signif((max(xdesign)-min(xdesign))/2,digits=3))
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      # to avoid 'yellow'
      cols<-cols[-7]
    }
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weigth Matrices",ylim=c(min(xdesign),max(xdesign)+(max(xdesign)-min(xdesign))/4),ylab="",xaxt="n",yaxt="n",pch=4)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1])        
    for(i in 2:nx)
    {  
      points(xaxis1,optmatrix[,i],col=cols[i],pch=4)
      lines(xaxis1,optmatrix[,i],col=cols[i])
    }
    legend(0.5,max(xdesign)+(max(xdesign)-min(xdesign))/4,names(data)[1:nx],col=cols[1:nx],bty="n",lwd=1)
  
  
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }
  
    # right Plot
  
    par(fig=c(0.5,0.95,0.15,0.85),new=TRUE,bty="l",pty="s",las=1)
    plot(xaxis1,reoptplot[,1],xlab="Weight Matrices",ylab="",ylim=c(0,1.25+max(devstand)),xaxt="n",yaxt="n",pch=4)
    mtext("Predicted Response",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=c(0.25,0.625,1),labels=c("","",""))
    lines(xaxis1,reoptplot[,1],col=cols[nx+1])
    for(i in 2:ny)
    {
      points(xaxis1,reoptplot[,i],col=cols[nx+i],pch=4)
      lines(xaxis1,reoptplot[,i],col=cols[nx+i])
    }
    for(i in 1:ny)
    {
      mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols[nx+i],cex=0.6,line=1.2)
    }
    #deviation
    for(i in 1:numbW)
    {
      point1<-NULL
      point2<-NULL
      point1<-c(i,i)
      for(j in 1:ny)
      {
        point2[1]<-reoptplot[i,j]+devstand[i,j]
        point2[2]<-reoptplot[i,j]-devstand[i,j]
        points(point1,point2,col=cols[nx+j],pch=4)
        lines(point1,point2,col=cols[nx+j])
      }
    }
    nam<-names(data)[(nx+1):(nx+ny)]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=") 
    }
    legend(numbW*(0.55),1.25+max(devstand),nam,col=cols[(nx+1):(nx+ny)],bty="n",lwd=1)#1.25+max(devstand)
    cat("...Plot is done!\n")
    dimnames(optmatrix)<-list(xaxis1names,names(data)[1:nx])
    dimnames(reoptmatrix)<-list(xaxis1names,names(data)[(nx+1):(nx+ny)])
    dimnames(deviation)<-list(xaxis1names,names(data)[(nx+1):(nx+ny)])
    optval<-as.matrix(optval)
    dimnames(optval)<-list(xaxis1names,c(""))
    optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau) 
    names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValue")
  return(optimres)
  }


