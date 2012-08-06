\name{locate}
\alias{locate}
\title{Selection of a compromise}
\description{
The function
\code{locate} allows the user to choose a point as a good compromise on the right plot and \code{locate} returns 
the corresponding design parameters. 
}
\usage{
locate(x,xlu=NaN,no.col=FALSE,standard=TRUE,col=1,lty=1,
bty="l",pty="s",las=1,adj=0.5,cex=1,cex.lab=0.8,cex.axis=0.8,
xlab=c("Stretch Vector","Stretch Vector"),ylab=c("Parameter Setting","Predicted Response"),...)
}
\arguments{
 \item{x}{object from JOP 
}
  \item{xlu}{The x-coordinate that indicates where the user assumes a good compromise,
  see Details  
  }
    \item{no.col}{ 
  If TRUE the plot will be gray scaled. Otherwise the plot will be coloured.
  }
    \item{standard}{ 
  If TRUE the standard deviations will be displayed on the right hand plot. 
  }
  \item{col}{Graphical argument, see details.
  }
    \item{lty}{Graphical argument, see details.
  }
    \item{xlab}{Graphical argument, see details. 
  }
    \item{ylab}{Graphical argument, see details.
  }
  \item{bty,pty,las,cex,adj,cex.lab,cex.axis}{Graphical arguments
  }
    \item{...}{ 
  Further graphical arguments passed to \code{\link{plot}}. 
  }
}

\value{
locate returns a list containing the following elements:
\item{ChosenResponses}{selected responses by user} 
\item{ChosenParameters}{corresponding selected parameters}  
}
\details{
The function \code{locate} asks the user to choose a compromise on the right hand plot 
and it displays the nearest calculated points by means of vertical lines. Furthermore it returns 
the chosen responses together with the corresponding parameters.

Let nx be the number of parameters (number of columns of datax) and ny be the number
of responses (number of columns of datay). Then col and lty must have length nx+ny. Otherwise
predefined grey colors (for no.col=TRUE) or standard colors 1, 2, ..., nx+ny are used. 
The arguments xlab and ylab must have length two, where the first entry contains the label
for x-axis and y-axis of the left hand plot and the second entry contains the label
for x-axis and y-axis of the right hand plot. Additional graphical arguments can be plugged in.
}
\references{
Sonja Kuhnt and Martina Erdbruegge (2004), A strategy of robust paramater design for multiple responses, 
Statistical Modelling; 4: 249-264, TU Dortmund.

Martina Erdbruegge, Sonja Kuhnt and Nikolaus Rudak (2011), Joint optimization of independent multiple responses based on loss functions,
Quality and Reliability Engineering International 27, doi: 10.1002/qre.1229.

Joseph J. Pignatiello (1993), Strategies for robust multiresponse quality engineering, IIE Transactions 25, 5-15, Texas A M University.

Alexios Ghalanos and Stefan Theussl (2010). Rsolnp: General Non-linear Optimization Using Augmented Lagrange Multiplier Method. R package version 1.0-6. 

Mathieu Ribatet and Bertrand Iooss (2010). JointModeling: Joint Modelling of Mean and Dispersion. R package version 1.0-2.
}
\author{
Sonja Kuhnt and Nikolaus Rudak
}


\examples{
# Example: Sheet metal hydroforming process
data(datax)
data(datay)
outtest<-JOP(numbW=5,tau=list(0,0.05),datax=datax,datay=datay)
locate(outtest,xlu=3)
}

