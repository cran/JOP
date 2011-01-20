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
  
  
       ####################################################################
       ####################################################################
       #########      Help Function to find the least      ################
       #########      significant parameter                ################
       ####################################################################
       ####################################################################
       
indmin<-function(x,y,wert)
{
    # (y)
    k<-0
    a<-wert
    ind<-0
    w<-NULL
    x<-abs(x)
    if(length(x)>1)
    {
      for(i in 2:length(x))
      {
        if(x[i]<=wert)
        {
          if(x[i]<a)
          {
            a<-x[i]
            ind<-i
            k<-1
          }
        }
      }
      if(k==1)
      {
        w<-y[-ind]
      }
    }
  return(list(a,ind,k,w))
}


       ####################################################################
       ####################################################################
       #########      Help Functions for creating formulas       ##########
       ####################################################################
       ####################################################################

       #################
       ## Mean-Modell ##
       #################

crformM<-function(resp,x)
{
   if(length(x)==1)
   {
    form<-as.formula(paste(resp,paste("~","1",sep=""),sep=""))
   }
   else
   {
    form<-as.formula(paste("'",resp,paste("~",paste(x[-1],collapse="+",sep=""),sep=""),"'",sep=""))
   }
   return(form)
}
      #################
      ## Disp-Modell ##
      #################

crformD<-function(x)
{
   if(length(x)==1)
   {
    form<-as.formula(paste("'d~","1'",sep=""))
   }
   else
   {
   form<-as.formula(paste("'d~",paste(x[-1],collapse="+",sep=""),"'",sep=""))
   return(form)
   }
}

        
       ###############################################
       ###############################################
       #########       Function 'JOP'       ##########
       ###############################################
       ###############################################
    
  JOP<-function(nx=2,ny=1,Wstart=1,Wend=1,numbW=1,d=c(1,0),optreg=0,tau=NULL,interact=0,
                quad=0,main.disp=0,interact.disp=0,quad.disp=0,drplot=TRUE,data=NULL,mean.model=NULL,var.model=NULL,solver=0,no.col=FALSE,standard=TRUE)
  {
  ## solver=0: nlminb
  ## solver=1: rgenoud
  ## optreg=0: Optimization region is a sphere
  ## optreg=1: Optimization region is a cube
  ## mean.model, var.model : List of functions of type mean.model<-function(x)
  ##                           that return a single response value
  # Input:
  # nx:                  Number of variables
  # ny:                  Number of responses
  # data:                DataFrame
  # Wstart,Wend,numbW,d: Parameter weight matrices
  #                      numbW: number of weight matrices
  #                      Formula: log wt= d*log at
  #                      Where: log at=seq(Wstart,Wend,(Wstart-Wend)/numbW)


  if(nx<2)
  {
    return("More then 1 variable required!!")
  }

      #############################################
      #############################################
      ##### Modeling by double gen lin models #####
      #############################################
      #############################################
      
  shifter<-0
  if(is.null(mean.model) & is.null(var.model))
  {
    shifter<-1
    print("Automatic Modeling starts...\n")
    print("\n")
        
    namen<-names(data)
  
    formeln<-paste("fm",1:ny,sep="")
  
    flist<-vector("list",ny)
    for(i in 1:ny)
    {
     names(flist)[i]<-formeln[i] 
    }
  
    xnames<-namen[1:nx]
  
    resp<-namen[(nx+1):(nx+ny)]
    if(interact==0)
    {
      rights<-paste(xnames,collapse="+",sep="")
    }
    if(interact==1)
    {
      rights<-paste("(",paste(xnames,collapse="+",sep=""),")^2",sep="")
    }
    if(interact!=0 && interact!=1)
    {
       return("Choose interact=1 or interact=0 please!!")
    }
    if(quad==1)
    {
      quadnames<-NULL
      for(i in 1:nx)
      {
        quadnames[i]<-paste("I(",xnames[i],"^2)",sep="")  
      }
      quadn<-paste(quadnames,collapse="+",sep="")
      rights<-paste(rights,quadn,sep="+")
    }
    
    if(interact.disp==0)
    {
      rights.disp<-paste(xnames,collapse="+",sep="")
    }
    if(interact.disp==1)
    {
      rights.disp<-paste("(",paste(xnames,collapse="+",sep=""),")^2",sep="")
    }
    if(interact.disp!=0 && interact.disp!=1)
    {
       return("Choose interact=1 or interact=0 please!!")
    }
    if(quad.disp==1)
    {
      quadnames.disp<-NULL
      for(i in 1:nx)
      {
        quadnames.disp[i]<-paste("I(",xnames[i],"^2)",sep="")  
      }
      quadn.disp<-paste(quadnames.disp,collapse="+",sep="")
      rights.disp<-paste(rights.disp,quadn.disp,sep="+")
    }
    if(interact.disp==0 && quad.disp==0 && main.disp==0)
    {
      rights.disp<-1
    }
    lefts<-paste(resp," ~ ")
  
        
        #################
        ## Mean-Modell ##
        #################
  
  
    for(i in 1:ny)
    {
      flist[[i]]<-paste("'",lefts[i],rights,"'",sep="")
    }

  
  
  
        #######################
        ## Dispersion-Modell ##
        #######################
    
    dispf<-vector("list",ny)
    for(i in 1:ny)
    {
      dispf[[i]]<-paste("'d~",rights.disp,"'",sep="")
    }

  
    cat("\n")
    cat("Automatic Modeling starts...\n")
    cat("\n")
    
    outmod<-vector("list",ny)
    namenlist<-as.list(namen[(nx+1):(nx+ny)])
    names(outmod)<-namenlist
    
    wert<-qt(0.95,dim(data)[1]-nx-1)
    for(i in 1:ny)
    {
      k<-1
      
      while(k==1)
      {
        daten<-data[c(nx+i,1:nx)]
        outmod[[i]]<-fitjoint("glm",flist[[i]],dispf[[i]],data=daten)
        outm<-summary(outmod[[i]]$mod.mean)
        outd<-summary(outmod[[i]]$mod.disp)
        
        helpa<-outm$coefficients[,3] 
        helpb<-outd$coefficients[,3] 

        minM<-indmin(helpa,names(helpa),wert)
        minD<-indmin(helpb,names(helpb),wert)
        
        if(minM[[3]]==0 & minD[[3]]==0)
        {
          k<-0 
        }
        else
        {
          if(minM[[3]]==0)
          {
              lxnamesd<-minD[[4]]
              dispf[[i]]<-crformD(lxnamesd)
          }
          if(minD[[3]]==0)
          {
              lxnames<-minM[[4]]
              flist[[i]]<-crformM(resp[i],lxnames)
          }
          if(minM[3]!=0 & minD[3]!=0)
          {
            if(minM[[1]]<minD[[1]])
            {
              lxnames<-minM[[4]]
              flist[[i]]<-crformM(resp[i],lxnames)
            }
            else
            {
              lxnamesd<-minD[[4]]
              dispf[[i]]<-crformD(lxnamesd)
            }
          }
        }
      }
      lxnames<-NULL
      lxnamesd<-NULL
    }
    # If the user only wants to calculate models, 
    # then outmod is returned
    if(drplot==FALSE)
    {
      return(outmod)
    }
    
    if(is.null(data))
    {
      return("data required!!")
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
      return("Choose a solver! '0' for 'nlminb' or '1' for goslnp or '2' for genoud....")
    } 
    ### Function for Mean and Disp
    mean.model<-vector("list",ny)
    var.model<-vector("list",ny)

    meanmd<-function(x,p)
    {
      allx<-NULL
      xnames<-names(data)[1:nx]
      coeff<-outmod[[p]]$mod.mean$coefficients
      value<-coeff[1]
      if(length(coeff)>1)
      {
        for(i in 2:length(coeff))
        {
           for(j in 1:length(xnames))  
           {
             if(names(coeff)[i]==xnames[j])
             {
                value<-value+coeff[i]*x[j]
             }
             if(names(coeff)[i]==paste("I(",xnames[j],"^2)",sep=""))
             {
                value<-value+coeff[i]*x[j]^2
             }
             for(k in j:length(xnames))
             {
                if(k!=j && names(coeff)[i]==paste(xnames[j],":",xnames[k],sep="") || coeff[i]==paste(xnames[k],":",xnames[j],sep=""))
                {
                  value<-value+coeff[i]*x[j]*x[k]
                }
             }
           }
        }
      }
      names(value)<-NULL
      return(value)
    }
    varmd<-function(x,p)
    {
      allx<-NULL
      xnames<-names(data)[1:nx]
      coeff<-outmod[[p]]$mod.disp$coefficients
      value<-coeff[1]
      if(length(coeff)>1)
      {
        for(i in 2:length(coeff))
        {
         for(j in 1:length(xnames))  
         {
           if(names(coeff)[i]==xnames[j])
           {
              value<-value+coeff[i]*x[j]
           }
           if(names(coeff)[i]==paste("I(",xnames[j],"^2)",sep=""))
           {
              value<-value+coeff[i]*x[j]^2
           }
           for(k in j:length(xnames))
           {
              if(k!=j && names(coeff)[i]==paste(xnames[j],":",xnames[k],sep="") || coeff[i]==paste(xnames[k],":",xnames[j],sep=""))
              {
                value<-value+coeff[i]*x[j]*x[k]
              }
           }
         }
        }
      }
      names(value)<-NULL
      return(exp(value))
    }
  
    cat("\n")
    cat("...Functions for mean and dispersion build...\n")
    cat("\n")         
  
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
         cat("...weight matrices build...\n")
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
               zw[j]<-as.numeric(varmd(xdesign[j,],i))  
           }
           diagA[i]<-1/sqrt(mean(zw))
           zw<-NULL 
         }  
         A<-diag(diagA)
         print(A)
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
        normvec<-NULL
        for(i in 1:dim(data)[1])
        {
          normvec[i]<-sqrt(t(as.numeric(data[i,1:nx]))%*%as.numeric(data[i,1:nx]))
        }
        upper<-c(max(normvec),rep(pi,nx-1))
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
               varval[j]<-as.numeric(varmd(trafopar(x),j))
               meanval[j]<-as.numeric(meanmd(trafopar(x),j))    
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
              reoptmatrix[i,k]<-as.numeric(meanmd(trafopar(opt$par),k))
              deviation[i,k]<-as.numeric(varmd(trafopar(opt$par),k))
            }
          }
          if(solver==1)
          { 
            opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2,control=list(trace=0))
            optval[i]<-opt$values[length(opt$values)]
            # The results are stored in optmatrix, reoptmatrix and deviation
            optmatrix[i,]<-trafopar(opt$pars)
            for(k in 1:ny)
            {
              reoptmatrix[i,k]<-as.numeric(meanmd(trafopar(opt$pars),k))
              deviation[i,k]<-as.numeric(varmd(trafopar(opt$pars),k))
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
              reoptmatrix[i,k]<-as.numeric(meanmd(trafopar(opt$par),k))
              deviation[i,k]<-as.numeric(varmd(trafopar(opt$par),k))
            }
          }
                                                       
        }
      }
      if(optreg==1)
      {
        lower<-NULL
        upper<-NULL
        for(i in 1:nx)
        {
          lower[i]<-min(data[,i])
          upper[i]<-max(data[,i])
        }
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
               varval[j]<-as.numeric(varmd(x,j))
               meanval[j]<-as.numeric(meanmd(x,j))    
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
              reoptmatrix[i,k]<-as.numeric(meanmd(opt$par,k))
              deviation[i,k]<-as.numeric(varmd(opt$par,k))
            }
          }
          if(solver==1)
          { 
            opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2,control=list(trace=0))
            optval[i]<-opt$values[length(opt$values)]
            # The results are stored in optmatrix, reoptmatrix and deviation
            optmatrix[i,]<-opt$pars
            for(k in 1:ny)
            {
              reoptmatrix[i,k]<-as.numeric(meanmd(opt$pars,k))
              deviation[i,k]<-as.numeric(varmd(opt$pars,k))
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
              reoptmatrix[i,k]<-as.numeric(meanmd(opt$par,k))
              deviation[i,k]<-as.numeric(varmd(opt$par,k))
            }
          }
                                                       
        }
      }
    }
    else
    {
      if(!is.list(mean.model) && !is.list(var.model))
      {
        return("mean.model and var.model have to be of type 'list'!")
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
        normvec<-NULL
        for(i in 1:dim(data)[1])
        {
          normvec[i]<-sqrt(t(as.numeric(data[i,1:nx]))%*%as.numeric(data[i,1:nx]))
        }
        upper<-c(max(normvec),rep(pi,nx-1))
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
            opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2,control=list(trace=0))
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
        lower<-NULL
        upper<-NULL
        for(i in 1:nx)
        {
          lower[i]<-min(data[,i])
          upper[i]<-max(data[,i])
        }
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
            opt<-gosolnp(fun=riscfun,LB=lower,UB=upper,n.restarts=2,control=list(trace=0))
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
    }
    deviation<-sqrt(deviation)
    cat("...Optimal Parameters and Responses received...\n")
    cat("\n")
    if(no.col==FALSE)
    {
    #################################
    #################################
    #### Joint optimization Plot ####
    #################################
    #################################
    
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      # to avoid 'yellow'
      cols<-cols[-7]
    }     
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weigth Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1])        
    for(i in 2:nx)
    {  
      points(xaxis1,optmatrix[,i],col=cols[i],pch=NA)
      lines(xaxis1,optmatrix[,i],col=cols[i])
    }
    legend(0.5,max(abs(optmatrix))+0.25*max(abs(optmatrix)),names(data)[1:nx],col=cols[1:nx],bty="n",lwd=1)
  
    
    
    # right Plot
    
      
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      if(max(reoptmatrix[,i])>tau[i])
      {
       reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
       devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
      }
      else
      {
       reoptplot[,i]<-reoptmatrix[,i]/tau[i]
       devstand[,i]<-deviation[,i]/tau[i] 
      }
    }
  

  
    par(fig=c(0.5,0.95,0.15,0.85),new=TRUE,bty="l",pty="s",las=1)
    plot(xaxis1,reoptplot[,1],xlab="Weight Matrices",ylab="",ylim=c(0,1.25+max(devstand)),xaxt="n",yaxt="n",pch=NA)
    mtext("Predicted Response",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=c(0.25,0.625,1),labels=c("","",""))

    

    #deviation
    if(standard==TRUE)
    {
      for(j in 1:ny)
      {
        point11<-NULL
        point12<-NULL
        point21<-NULL
        point22<-NULL
        for(i in 1:numbW)
        {
          point11[i]<-i
          point21[i]<-numbW-i+1
          point22[i]<-reoptplot[numbW-i+1,j]+devstand[numbW-i+1,j]
          point12[i]<-reoptplot[i,j]-devstand[i,j]
          #points(point1,point2,col=cols[nx+j],pch=NA)
          #lines(point1,point2,col=cols[nx+j])
        }
        point1<-c(point11,point21)
        point2<-c(point12,point22)
        polygon(point1,point2,col=cols[nx+j],border=FALSE)
      }    
    
    lines(xaxis1,reoptplot[,1],col="black")#cols[nx+1])                            
    for(i in 2:ny)
    {
      points(xaxis1,reoptplot[,i],pch=NA)
      lines(xaxis1,reoptplot[,i],col="black")#cols[nx+i])
    }
    for(i in 1:ny)
    {
      mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols[nx+i],cex=0.6,line=1.2)
    }
    }
    if(standard==FALSE)
    {
      lines(xaxis1,reoptplot[,1],col=cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col=cols[nx+i])
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols[nx+i],cex=0.6,line=1.2)
      }
    }
    ######################
    ### Target Values ####
    ######################
        
    axis4<-NULL
    for(i in 1:ny)
    {
      axis4[i]<-tau[i]/max(reoptmatrix[,i])
    }
    zaehler<-vector("list",length(tau))
    for(i in 1:length(tau))
    {
    counter<-NULL
      for(j in 1:length(tau))
      {
        if(abs(axis4[i]-axis4[j])<0.065 && i!=j)
        {
          counter[i]<-j
        }
      }
      zaehler[[i]]<-sort(c(i,counter))
    }
    for(i in 1:length(tau))
    {
      point1<-NULL
      point2<-NULL
      for(j in 1:length(zaehler[[i]]))
      {
        if(i==zaehler[[i]][j])
        {
          point1<-c(1+(j-1)*(numbW-1)/length(zaehler[[i]]),1+j*(numbW-1)/length(zaehler[[i]]))
          point2<-c(axis4[i],axis4[i])
          lines(point1,point2,col=cols[nx+i])
          mtext(signif(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols[nx+i],cex=0.6,line=1.2)
        }
      }
    }  

    ##################################
    ####### End: Target Values #######
    ##################################
    
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
    if(shifter==1)
    {
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau,outmod)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValue","DGLM") 
    }
    else
    {
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValue") 
    }
    
    }
    else
    {
    #################################
    #################################
    #### Joint optimization Plot ####
    #################################
    #################################
    
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
    cols1<-gray(seq(0.5,0.7,length=nx))
    cols2<-gray(seq(0.5,0.7,length=ny))
    cols3<-gray(seq(0.75,0.95,length=ny))
    
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weigth Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1])        
    for(i in 2:nx)
    {  
      points(xaxis1,optmatrix[,i],pch=NA)
      lines(xaxis1,optmatrix[,i],col=cols1[i])
    }
    legend(0.5,max(abs(optmatrix))+0.25*max(abs(optmatrix)),names(data)[1:nx],col=cols1[1:nx],bty="n",lwd=1)
  
    
    
    # right Plot
    
      
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      if(max(reoptmatrix[,i])>tau[i])
      {
       reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
       devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
      }
      else
      {
       reoptplot[,i]<-reoptmatrix[,i]/tau[i]
       devstand[,i]<-deviation[,i]/tau[i] 
      }
    }
  

  
    par(fig=c(0.5,0.95,0.15,0.85),new=TRUE,bty="l",pty="s",las=1)
    plot(xaxis1,reoptplot[,1],xlab="Weight Matrices",ylab="",ylim=c(0,1.25+max(devstand)),xaxt="n",yaxt="n",pch=NA)
    mtext("Predicted Response",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=c(0.25,0.625,1),labels=c("","",""))

    

    #deviation
    if(standard==TRUE)
    {
      for(j in 1:ny)
      {
        point11<-NULL
        point12<-NULL
        point21<-NULL
        point22<-NULL
        for(i in 1:numbW)
        {
          point11[i]<-i
          point21[i]<-numbW-i+1
          point22[i]<-reoptplot[numbW-i+1,j]+devstand[numbW-i+1,j]
          point12[i]<-reoptplot[i,j]-devstand[i,j]
          #points(point1,point2,col=cols[nx+j],pch=NA)
          #lines(point1,point2,col=cols[nx+j])
        }
        point1<-c(point11,point21)
        point2<-c(point12,point22)
        polygon(point1,point2,col=cols3[j],border=FALSE)
      }    

    
      lines(xaxis1,reoptplot[,1],col="black")#cols2[1])#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col="black")#cols2[i])
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols2[i],cex=0.6,line=1.2)
      }
    }
    if(standard==FALSE)
    {
      lines(xaxis1,reoptplot[,1],col=cols2[1])#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col=cols2[i])
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols2[i],cex=0.6,line=1.2)
      }
    }
    ######################
    ### Target Values ####
    ######################
        
    axis4<-NULL
    for(i in 1:ny)
    {
      axis4[i]<-tau[i]/max(reoptmatrix[,i])
    }
    zaehler<-vector("list",length(tau))
    for(i in 1:length(tau))
    {
    counter<-NULL
      for(j in 1:length(tau))
      {
        if(abs(axis4[i]-axis4[j])<0.065 && i!=j)
        {
          counter[i]<-j
        }
      }
      zaehler[[i]]<-sort(c(i,counter))
    }
    for(i in 1:length(tau))
    {
      point1<-NULL
      point2<-NULL
      for(j in 1:length(zaehler[[i]]))
      {
        if(i==zaehler[[i]][j])
        {
          point1<-c(1+(j-1)*(numbW-1)/length(zaehler[[i]]),1+j*(numbW-1)/length(zaehler[[i]]))
          point2<-c(axis4[i],axis4[i])
          lines(point1,point2,col=cols2[i])
          mtext(signif(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols2[i],cex=0.6,line=1.2)
        }
      }
    }  

    ##################################
    ####### End: Target Values #######
    ##################################
    
    nam<-names(data)[(nx+1):(nx+ny)]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=") 
    }
    legend(numbW*(0.55),1.25+max(devstand),nam,col=cols2,bty="n",lwd=1)#1.25+max(devstand)
    cat("...Plot is done!\n")
    dimnames(optmatrix)<-list(xaxis1names,names(data)[1:nx])
    dimnames(reoptmatrix)<-list(xaxis1names,names(data)[(nx+1):(nx+ny)])
    dimnames(deviation)<-list(xaxis1names,names(data)[(nx+1):(nx+ny)])
    optval<-as.matrix(optval)
    dimnames(optval)<-list(xaxis1names,c(""))
    if(shifter==1)
    {
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau,outmod)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValue","DGLM") 
    }
    else
    {
      optimres<-list(optmatrix,reoptmatrix,deviation,optval,tau)
      names(optimres)<-list("Parameters","Responses","StandardDeviation","OptimalValue","TargetValue") 
    }
    }
    return(optimres)
}
      
      