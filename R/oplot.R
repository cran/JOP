.packageName<-'JOP'

oplot <-
function(out,no.col=FALSE,standard=TRUE)
{  
  Wstart<-out$ValW[1]
  Wend<-out$ValW[2]
  ##  out is output of 'JOP'
  if(no.col==FALSE)
  {
    ## Warning Messages:
    if(is.null(out))
    {
      return("Run 'JOP' first! 'oplot' needs output of 'JOP'!")
    }
  
    ## Setting Values
    nx<-dim(out[[1]])[2]
    ny<-dim(out[[2]])[2]
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      cols<-cols[-7]
    }
    numbW<-dim(out[[2]])[1]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    yaxis1<-seq(-round(max(abs(optmatrix)),digits=3),round(max(abs(optmatrix)),digits=3),round(max(abs(optmatrix)),digits=3))
  
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Stretch Vector",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend))
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1])
    for(i in 2:nx)
    {
      points(xaxis1,optmatrix[,i],col=cols[i],pch=NA)
      lines(xaxis1,optmatrix[,i],col=cols[i])
    }
    legend("topright",dimnames(out$Parameters)[[2]][1:nx],col=cols[1:nx],bty="n",lwd=1)
  
  
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }
  
    # right Plot
  
    par(fig=c(0.5,0.95,0.15,0.85),new=TRUE,bty="l",pty="s",las=1)
    plot(xaxis1,reoptplot[,1],xlab="Stretch Vector",ylab="",ylim=c(0,1.25+max(devstand)),xaxt="n",yaxt="n",pch=NA)
    mtext("Predicted Response",side=3,at=1,cex=0.6)
    axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend))
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
        mtext(c(round(0.25*max(reoptmatrix[,i]),digits=2),round(0.625*max(reoptmatrix[,i]),digits=2),round(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols[nx+i],cex=0.6,line=1.2)
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
        mtext(c(round(0.25*max(reoptmatrix[,i]),digits=2),round(0.625*max(reoptmatrix[,i]),digits=2),round(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols[nx+i],cex=0.6,line=1.2)
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
            mtext(round(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols[nx+i],cex=0.6,line=1.2)
          }
        }
      }  
  
      ##################################
      ####### End: Target Values #######
      ##################################
      
    nam<-dimnames(out$Responses)[[2]]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=")
    }
    legend("topright",nam,col=cols[(nx+1):(nx+ny)],bty="n",lwd=1)
  }
  else
  {
      ## Warning Messages:
    if(is.null(out))
    {
      return("Run 'JOP' first! 'oplot' needs output of 'JOP'!")
    }
  
    ## Setting Values
    nx<-dim(out[[1]])[2]
    ny<-dim(out[[2]])[2]
    cols<-1:(nx+ny+1)
    if((nx+ny)>=7)
    {
      cols<-cols[-7]
    }
    numbW<-dim(out[[2]])[1]
    optmatrix<-out[[1]]
    reoptmatrix<-out[[2]]
    deviation<-out[[3]]
    tau<-out[[5]]
  
    ## Values to label the axes
    xaxis1<-1:numbW
    xaxis2<-1:(numbW+ny-1)
    xaxis1names<-paste("W",xaxis1,sep="")
    yaxis1<-seq(-round(max(abs(optmatrix)),digits=3),round(max(abs(optmatrix)),digits=3),round(max(abs(optmatrix)),digits=3))
    cols1<-gray(seq(0.5,0.6,length=nx))
    cols3<-gray(seq(0.5,0.6,length=ny))
    cols2<-gray(seq(0.85,0.95,length=ny))
    
    # left plot
    par(fig=c(0,0.45,0.15,0.85),lwd=1,lty=6,bty="l",pty="s",las=1,cex=0.6,adj=0.5)
    plot(xaxis1,optmatrix[,1],xlab="Stretch Vector",ylim=c(-max(abs(optmatrix)),max(abs(optmatrix))+0.25*max(abs(optmatrix))),ylab="",xaxt="n",yaxt="n",pch=NA)
    mtext("Parameter Setting",side=3,at=1,cex=0.6)
    axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend))
    axis(2,at=yaxis1)
    lines(xaxis1,optmatrix[,1],lty=1,col=cols3[1])        
    for(i in 2:nx)
    {  
      points(xaxis1,optmatrix[,i],pch=NA)
      lines(xaxis1,optmatrix[,i],lty=i,col=cols3[i])
    }
    legend("topright",dimnames(out$Parameters)[[2]],lty=1:nx,bty="n",lwd=1)
  
    
    
    # right Plot
    
    reoptplot<-matrix(NaN,ncol=ny,nrow=numbW)
    devstand<-matrix(NaN,ncol=ny,nrow=numbW)
    for(i in 1:ny)
    {
      reoptplot[,i]<-reoptmatrix[,i]/max(reoptmatrix[,i])
      devstand[,i]<-deviation[,i]/max(reoptmatrix[,i])
    }

  

  
    par(fig=c(0.5,0.95,0.15,0.85),new=TRUE,bty="l",pty="s",las=1)
    plot(xaxis1,reoptplot[,1],xlab="Stretch Vector",ylab="",ylim=c(0,1.25+max(devstand)),xaxt="n",yaxt="n",pch=NA)
    mtext("Predicted Response",side=3,at=1,cex=0.6)
    axis(1,at=c(1,1+(numbW-1)/2,numbW),labels=c(Wstart,Wstart+0.5*(Wend-Wstart),Wend))
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
          
      lines(xaxis1,reoptplot[,1],lty=1,col=cols3[1])#cols2[1])#cols[nx+1])                            
      for(i in 2:ny)
      {
        points(xaxis1,reoptplot[,i],pch=NA)
        lines(xaxis1,reoptplot[,i],lty=i,col=cols3[i])#=cols2[i])
      }
      for(i in 1:ny)
      {
        mtext(c(round(0.25*max(reoptmatrix[,i]),digits=2),round(0.625*max(reoptmatrix[,i]),digits=2),round(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
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
        mtext(c(round(0.25*max(reoptmatrix[,i]),digits=2),round(0.625*max(reoptmatrix[,i]),digits=2),round(max(reoptmatrix[,i]),digits=2)),side=2,at=c(0.25-(ny-1)*0.035+(i-1)*0.065,0.625-(ny-1)*0.035+(i-1)*0.065,1-(ny-1)*0.045+(i-1)*0.065),col=cols3[i],cex=0.6,line=1.2)
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
          mtext(round(tau[i],digits=2),side=4,at=axis4[i]-(j-1)*0.065,col=cols3[i],cex=0.6,line=1.2)
        }
      }
    }  

    ##################################
    ####### End: Target Values #######
    ##################################
    
    nam<-dimnames(out$Responses)[[2]]
    for(i in 1:length(nam))
    {
      nam[i]<-paste(paste(paste(nam[i],";",sep=""),"target",sep=" "),tau[i],sep="=") 
    }
    legend("topright",nam,col=cols3,lty=1:ny,bty="n",lwd=1)#1.25+max(devstand)
  }
  return("Plot is done!")
}


