\name{plot.JOP}
\alias{plot.JOP}
\title{Displaying the Joint Optimization Plot}
\description{
The function \code{plot.JOP} takes the output 
produced by \code{JOP} and returns the joint optimization plot.
}
\usage{
\method{plot}{JOP}(x, no.col = FALSE, standard = TRUE, col = 1, lty = 1, bty = "l",
las = 1 ,adj = 0.5 ,cex = 1 ,cex.lab = 1 ,cex.axis = 1,
xlab = c("Stretch Vector" , "Stretch Vector"),
ylab = c("Parameter Setting" , "Predicted Response"),lwd=1,...)
}
\arguments{
 \item{x}{object from JOP 
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
  \item{bty,las,cex,adj,cex.lab,cex.axis,lwd}{Graphical arguments
  }
    \item{...}{ 
  Further graphical arguments passed to \code{\link{plot}}. 
  }
}
\details{
Let nx be the number of parameters (number of columns of datax) and ny be the number
of responses (number of columns of datay). Then col and lty must have length nx+ny. Otherwise
predefined grey colors (for no.col=TRUE) or standard colors 1, 2, ..., nx+ny are used. 
The arguments xlab and ylab must have length two, where the first entry contains the label
for x-axis and y-axis of the left hand plot and the second entry contains the label
for x-axis and y-axis of the right hand plot. Additional graphical arguments can be plugged in.  
}
\references{
Sonja Kuhnt and Martina Erdbruegge (2004). A strategy of robust paramater design for multiple responses, 
Statistical Modelling; 4: 249-264, TU Dortmund.

Martina Erdbruegge, Sonja Kuhnt and Nikolaus Rudak (2011). Joint optimization of independent multiple responses based on loss functions,
Quality and Reliability Engineering International 27, doi: 10.1002/qre.1229.

Joseph J. Pignatiello (1993). Strategies for robust multiresponse quality engineering, IIE Transactions 25, 5-15, Texas A M University.

Alexios Ghalanos and Stefan Theussl (2012). Rsolnp: General Non-linear Optimization Using Augmented Lagrange Multiplier Method. R package version 1.12. 

Peter K Dunn and Gordon K Smyth (2012). dglm: Double generalized linear models, R package version 1.6.2.

Sonja Kuhnt, Nikolaus Rudak (2013). Simultaneous Optimization of Multiple Responses with the R Package JOP, 
Journal of Statistical Software, 54(9), 1-23, URL http://www.jstatsoft.org/v54/i09/.
}


\examples{
# Example: Sheet metal hydroforming process
outtest <- JOP(datax = datax, datay = datay, tau = list(0 , 0.05), numbW = 5)

# Several graphical parameters can be plugged in
plot(outtest, col = 5:8)
}

