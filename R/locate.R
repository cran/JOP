.packageName<-'locate'

locate <-
function(data,out,xlu=NaN,no.col=FALSE,standard=TRUE)
{
  numbW<-dim(out$Responses)[1]
  ##  out is output of 'jointplot'
  if(no.col==FALSE)
  {
    ## Warning Messages:
    if(is.null(data))
    {
      return("data set is required!")
    }
    if(is.null(out))
    {
      return("Run 'jointplot' first! 'oplot' needs output of 'jointplot'!")
    }
    if(!is.nan(xlu))
    {
      if(xlu<1 || xlu>numbW)
      {
        cat("Choose a x-coordinate between 1 and numbW")
        cat("\n")
        return("Call locate again!")
      }
    }
    else
    {
      cat("Choose your preferred point on the right plot please!\n")
      cat("\n")
      flush.console()
    }
    ## Setting Values
    nx<-dim(out$Parameters)[2]
    ny<-dim(out$Responses)[2]
      cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      cols<-cols[-7]
    }

    xdesign<-data[,1:nx]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
  
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weight Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
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
  
  
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }
  
    # right Plot
  
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
        #points(xaxis1,reoptplot[,i],col=cols[nx+i],pch=NA)
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
        #points(xaxis1,reoptplot[,i],col=cols[nx+i],pch=NA)
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
     
    legend(numbW*(0.55),1.25+max(devstand),nam,col=cols[(nx+1):(nx+ny)],bty="n",lwd=1)
  
  
    ####
    ####  Location of points
    ####
  
  
  
   
    if(is.nan(xlu))
    {
      locp<-locator(n=1)
      xloc1<-c(locp$x,locp$x)  
    }
    else
    {
      xloc1<-c(xlu,xlu)  
    }
      
    ## Chosen Points:
    xl<-xloc1[1]
    for(i in 1:(numbW-1))
    {
      if((i<=xl) && (xl<(i+1)))
      {
        optp<-(optmatrix[i+1,]-optmatrix[i,])*(xl-floor(xl))+optmatrix[i,]
        reoptp<-(reoptmatrix[i+1,]-reoptmatrix[i,])*(xl-floor(xl))+reoptmatrix[i,]    
      }
    }
    if(xl==numbW)
    {
      optp<-out$Parameters[numbW,]
      reoptp<-out$Responses[numbW,]
    }  
      ##  out is output of 'jointplot'
  
    ## Warning Messages:
    if(is.null(data))
    {
      return("data set is required!")
    }
    if(is.null(out))
    {
      return("Run 'jointplot' first! 'oplot' needs output of 'jointplot'!")
    }
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      cols<-cols[-7]
    }
    ## Setting Values
    nx<-dim(out$Parameters)[2]
    ny<-dim(out$Responses)[2]
    numbW<-dim(out$Responses)[1]
    xdesign<-data[,1:nx]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
  
    # left plot
    par(fig=c(0,0.45,0.15,0.85),new=TRUE,lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weight Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1])
    for(i in 2:nx)
    {
      points(xaxis1,optmatrix[,i],col=cols[i],pch=NA)
      lines(xaxis1,optmatrix[,i],col=cols[i])
    }
  
  
    yloc1<-c(-max(abs(optmatrix)),max(abs(optmatrix)))
    points(xloc1,yloc1,pch=NaN)
    lines(xloc1,yloc1,lwd=2,lty="solid")
    xlp<-rep(xloc1,nx)
    xlr<-rep(xloc1,ny)
    for(i in 1:nx)
    {
      points(xlp[i],optp[i],col=cols[i],cex=3)
    }
    legend(0.5,max(abs(optmatrix))+0.25*max(abs(optmatrix)),names(data)[1:nx],col=cols[1:nx],bty="n",lwd=1)
  
  
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }
  
  
    # right Plot
  
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
        #points(xaxis1,reoptplot[,i],col=cols[nx+i],pch=NA)
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
        #points(xaxis1,reoptplot[,i],col=cols[nx+i],pch=NA)
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
      
    yloc1<-c(0,1.25)
    points(xloc1,yloc1,pch=NaN)
    lines(xloc1,yloc1,lwd=2,lty="solid")
    nam<-names(data)[(nx+1):(nx+ny)]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=")
    }
    for(i in 1:ny)
    {
      points(xlp[i],reoptp[i]/max(reoptmatrix[,i]),col="black",cex=3)#cols[nx+i],cex=3)
    }
    legend(numbW*(0.55),1.25+max(devstand),nam,col=cols[(nx+1):(nx+ny)],bty="n",lwd=1)
  
    opt<-list(optp,reoptp)
  
    names(opt)<-list("ChosenParameters","ChosenResponses")

    return(opt)
  }
  else
  {
      ## Warning Messages:
    if(is.null(data))
    {
      return("data set is required!")
    }
    if(is.null(out))
    {
      return("Run 'jointplot' first! 'oplot' needs output of 'jointplot'!")
    }
    if(!is.nan(xlu))
    {
      if(xlu<1 || xlu>numbW)
      {
        cat("Choose a x-coordinate between 1 and numbW")
        cat("\n")
        return("Call locate again!")
      }
    }
    else
    {
      cat("Choose your preferred point on the right plot please!\n")
      cat("\n")
      flush.console()
    }
    ## Setting Values
    ## Setting Values
    nx<-dim(out[[1]])[2]
    ny<-dim(out[[2]])[2]
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      cols<-cols[-7]
    }
    numbW<-dim(out[[2]])[1]
    xdesign<-data[,1:nx]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
    cols1<-gray(seq(0.5,0.6,length=nx))
    cols3<-gray(seq(0.5,0.6,length=ny))
    cols2<-gray(seq(0.85,0.95,length=ny))
    
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weight Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1],lty=1,col=cols3[1])        
    for(i in 2:nx)
    {  
      points(xaxis1,optmatrix[,i],pch=NA)
      lines(xaxis1,optmatrix[,i],lty=i,col=cols3[i])
    }
    legend(0.5,max(abs(optmatrix))+0.25*max(abs(optmatrix)),names(data)[1:nx],lty=1:nx,bty="n",lwd=1)
  
    
    
    # right Plot
    
      
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
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
      polygon(point1,point2,col=cols2[j],border=FALSE)
      }
          
      lines(xaxis1,reoptplot[,1],lty=1,col=cols3[1])#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],lty=i,col=cols3[i])
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
      }
    }
    if(standard==FALSE)
    {
      lines(xaxis1,reoptplot[,1],col=cols3[1],lty=1)#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col=cols3[i],lty=i)
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
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
          lines(point1,point2,col=cols3[i],lty=i)
          mtext(signif(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols3[i],cex=0.6,line=1.2)
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
    legend(numbW*(0.55),1.25+max(devstand),nam,col=cols3,lty=1:ny,bty="n",lwd=1)#1.25+max(devstand)  
  
    ####
    ####  Location of points
    ####
  
  
  
   
    if(is.nan(xlu))
    {
      locp<-locator(n=1)
      xloc1<-c(locp$x,locp$x)  
    }
    else
    {
      xloc1<-c(xlu,xlu)  
    }
      
    ## Chosen Points:
    xl<-xloc1[1]
    for(i in 1:(numbW-1))
    {
      if((i<=xl) && (xl<(i+1)))
      {
        optp<-(optmatrix[i+1,]-optmatrix[i,])*(xl-floor(xl))+optmatrix[i,]
        reoptp<-(reoptmatrix[i+1,]-reoptmatrix[i,])*(xl-floor(xl))+reoptmatrix[i,]    
      }
    }
    if(xl==numbW)
    {
      optp<-out$Parameters[numbW,]
      reoptp<-out$Responses[numbW,]
    }  
      ##  out is output of 'jointplot'
  
    ## Warning Messages:
    if(is.null(data))
    {
      return("data set is required!")
    }
    if(is.null(out))
    {
      return("Run 'jointplot' first! 'oplot' needs output of 'jointplot'!")
    }
    ## Setting Values
    nx<-dim(out$Parameters)[2]
    ny<-dim(out$Responses)[2]
    numbW<-dim(out$Responses)[1]
    xdesign<-data[,1:nx]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    parameternames<-paste("X",1:nx,sep="")
    yaxis1<-seq(-signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3),signif(max(abs(optmatrix)),digits=3))
  
    # left plot
    par(fig=c(0,0.45,0.15,0.85),new=TRUE,lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Weight Matrices",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=xaxis1,labels=xaxis1names)
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1],lty=1,col=cols3[1])
    for(i in 2:nx)
    {
      points(xaxis1,optmatrix[,i],pch=NA)
      lines(xaxis1,optmatrix[,i],lty=i,col=cols3[i])
    }
  
  
    yloc1<-c(-max(abs(optmatrix)),max(abs(optmatrix)))
    points(xloc1,yloc1,pch=NaN)
    lines(xloc1,yloc1,lwd=2,lty=1)
    xlp<-rep(xloc1,nx)
    xlr<-rep(xloc1,ny)
    for(i in 1:nx)
    {
      points(xlp[i],optp[i],col="black",cex=3)
    }
    legend(0.5,max(abs(optmatrix))+0.25*max(abs(optmatrix)),names(data)[1:nx],lty=1:nx,col=cols3,bty="n",lwd=1)
  
  
    # right Plot
    
      
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }

  
    # right Plot
  
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
      polygon(point1,point2,col=cols2[j],border=FALSE)
      }
          
      lines(xaxis1,reoptplot[,1],col=cols3[1],lty=1)#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col=cols3[i],lty=i)
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
      }
    }
    if(standard==FALSE)
    {
      lines(xaxis1,reoptplot[,1],col=cols3[1],lty=1)#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],col=cols3[i],lty=i)
      }
      for(i in 1:ny)
      {
        mtext(c(signif(0.25*max(reoptmatrix[,i]),digits=2),signif(0.625*max(reoptmatrix[,i]),digits=2),signif(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
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
            lines(point1,point2,col=cols3[i],lty=i)
            mtext(signif(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols3[i],cex=0.6,line=1.2)
          }
        }
      }  
  
      ##################################
      ####### End: Target Values #######
      ##################################
    
      ##################################
      ####### End: Target Values #######
      ##################################
      
    yloc1<-c(0,1.25)
    points(xloc1,yloc1,pch=NaN)
    lines(xloc1,yloc1,lwd=2,lty="solid")
    nam<-names(data)[(nx+1):(nx+ny)]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=")
    }
    for(i in 1:ny)
    {
      points(xlp[i],reoptp[i]/max(reoptmatrix[,i]),col="black",cex=3)#cols[nx+i],cex=3)
    }
    legend(numbW*(0.55),1.25+max(devstand),nam,lty=1:ny,col=cols3,bty="n",lwd=1)
  
    opt<-list(optp,reoptp)
  
    names(opt)<-list("ChosenParameters","ChosenResponses")

    return(opt)
    }
}


