.packageName<-'JOP'

oplot <-
function(data,out)
{
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
  yaxis1<-seq(signif(min(xdesign),digits=3),signif(max(xdesign),digits=3),signif((max(xdesign)-min(xdesign))/2,digits=3))

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
  legend(numbW*(0.55),1.25+max(devstand),nam,col=cols[(nx+1):(nx+ny)],bty="n",lwd=1)
  return("Plot is done!")
}


