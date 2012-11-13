\name{JOP}
\alias{JOP}
\title{Main function to minimize the risc function of a sequence of cost matrices}
\description{
\code{JOP} calculates optimal design parameters associated with a given sequence of cost matrices
based on the minimization of a risk function introduced by Pignatiello (1993). Furthermore \code{JOP} 
visualizes the optimal design parameters and the appropriate predicted responses using the joint 
optimization plot introduced by Kuhnt and Erdbruegge (2004). 
}
\usage{
JOP(Wstart=-5,Wend=5,numbW=10,d=NULL,optreg="sphere",Domain=NULL,tau,datax,datay,form.mean=NULL,form.disp=NULL,
family.mean=gaussian(),dlink="log",mean.model=NULL,var.model=NULL,joplot=F,solver="solnp")
}
\arguments{
  \item{Wstart}{
  value to calculate the sequence of weight matrices (see Details)
}
  \item{Wend}{
  value to calculate the sequence of weight matrices (see Details)
  }
  \item{numbW}{
  value to calculate the sequence of weight matrices (see Details)
}
  \item{d}{
  a vector with values to calculate the sequence of weight matrices (see Details)
}
  \item{optreg}{
  User can choose the Optimization region. 
  
  \code{optreg="box"}: box constraints  
  
  \code{optreg="sphere"}: sphere 
}
  \item{Domain}{
  box constraints. Column 1 for lower contraints and Column 2 for upper contraints.
  Row i corresponds to Parameter i.
}
  \item{tau}{
  list of numeric target values for the corresponding responses or specification "min"
  for minimization or "max" for maximization
}
  \item{datax}{
  data set with parameter settings from an experimental design (data.frame). 
  Columns have to be named.
}
  \item{datay}{
  data set with responses resulting from an experimental design (data.frame).
  Columns have to be named.
}
  \item{form.mean}{
  list of formulas for mean of each response
}
  \item{form.disp}{
  list of formulas for dispersion of each response
}
  \item{family.mean}{
  \code{\link{family}} for the mean
}
  \item{dlink}{
  \code{list of names of link functions for the dispersion models}
  }
  \item{mean.model}{
  list of functions that model the mean for the corresponding response
}
  \item{var.model}{
  list of functions that model the variance for the corresponding response
}
  \item{joplot}{ 
  logical, if TRUE then the joint optimization plot is displayed. 
  }
  \item{solver}{ 
  Default is "solnp" for three different starting points. Alternatively,
  "gosolnp" is especially recommended for complex programs.  
  }
}

\value{
JOP returns a list containing the following elements:
\item{Parameters}{The i-th row of this matrix contains the optimal Parameter setting appropriate to the i-th weight matrix} 
\item{Responses}{The i-th row of this matrix contains the predicted Responses appropriate to the i-th weight matrix}  
\item{StandardDeviation}{The i-th row of this matrix contains the standard deviation value for each response}
\item{OptimalValue}{This vector contains the optimal value of the risk function for each optimal parameter setting}
\item{TargetValueJOP}{Contains the target values for the correspoding responses used internally by JOP}
\item{TargetValueUSER}{Contains the target values for the correspoding responses specified by the user}
\item{DGLM}{If no models assigned then the list DGLM contains the calculated models for the mean and dispersion for every response}
\item{RiskminimalParameters}{Parameters that minimize the squared sum of single risks among all calculated Parameters}
\item{RiskminimalResponses}{Responses for the risk minimal parameters}
\item{diagA}{Diagonal values of the standardization matrix}
\item{ValW}{Values for Wstart and Wend}
}
\details{
The main function \code{JOP} is a package for multiresponse optimization which aims to minimize a risc function 
for a prespecified sequence of cost matrices. This sequence of cost matrices is specified by the arguments 
Wstart, Wend, numbW and d. The user can plug in target values for the responses or set to the target value
to "min" or "max" in order to minimize or maximize the corresponding response. 

JOP needs models for the mean and dispersion of each response which can be 
plugged in by means three different possibilities. \itemize{\item{First, the user can pass the models for mean and dispersion
as lists of functions in the parameter vector through the arguments var.model and mean.model. } \item{Secondly, the user
can plug in a list of formulas for each response for the mean and dispersion via the arguments form.mean 
and form.disp. Furthermore, a suitable link
and distribution assumption can be specified both for the mean and dispersion} \item{Finally, if the user does not plug in neither formulas nor models then JOP calculates automatically
double generalized linear models by means of the function \code{\link{dglm}} from package dglm. 
Furthermore, JOP performs a backward selection, starting from the full model with main effects, interactions
and quadratic terms, and afterwards dropping the least significant covariate in each step.}} The data sets \code{datax} and \code{datay}
are needed for model building. Both \code{datax} and \code{datay} have to be data frames where \code{datax} contains an experimental design with
settings for each parameter columnwise and \code{datay} contains the experimental results columnwise for every response. Additionally,
the columns of the data sets should be named, as exemplary demonstrated by data(datax) and data(datay). The optimization is performed
by the procedure \code{\link{solnp}} out of the package Rsolnp. JOP returns an object of class "JOP" which can be visualized by means of
\code{\link{plot.JOP}}. Details on the JOP method can be found in Erdbruegge et al. (2011).
}
\references{
Sonja Kuhnt and Martina Erdbruegge (2004). A strategy of robust paramater design for multiple responses, 
Statistical Modelling; 4: 249-264, TU Dortmund.

Martina Erdbruegge, Sonja Kuhnt and Nikolaus Rudak (2011). Joint optimization of independent multiple responses based on loss functions,
Quality and Reliability Engineering International 27, doi: 10.1002/qre.1229.

Joseph J. Pignatiello (1993). Strategies for robust multiresponse quality engineering, IIE Transactions 25, 5-15, Texas A M University.

Alexios Ghalanos and Stefan Theussl (2012). Rsolnp: General Non-linear Optimization Using Augmented Lagrange Multiplier Method. R package version 1.12. 

Peter K Dunn and Gordon K Smyth (2012). dglm: Double generalized linear models, R package version 1.6.2.
}

\examples{
# Example: Sheet metal hydroforming process
data(datax)
data(datay)

# Run JOP without Model specification
outtest<-JOP(numbW=5,tau=list(0,0.05),datax=datax,datay=datay,joplot=TRUE)
}

