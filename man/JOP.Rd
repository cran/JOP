\name{JOP}
\alias{JOP}
\alias{oplot}
\alias{JOPexample}
\alias{JOPexampledata}
\alias{locate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
JOP, oplot and locate
}
\description{
\code{JOP} calculates optimal design parameters associated with a given sequence of weight matrices
based on the minimization of a risk function introduced by Pignatiello (1993). Furthermore \code{JOP} 
visualizes the optimal design parameters and the appropriate predicted responses using the joint 
optimization plot introduced by Kuhnt and Erdbruegge (2004). The function \code{oplot} takes the output 
produced by \code{JOP} and conducts the visualization without optimization. Furthermore the function
\code{locate} allows the user to choose a point as a "good" compromise on the right plot and \code{locate} returns 
the corresponding design parameters. 
}
\usage{
JOP(nx=2,ny=1,Wstart=-5,Wend=5,numbW=20,d=c(1,0),optreg=0,tau=NULL,interact=1,
quad=1,main.disp=0,interact.disp=0,quad.disp=0,data=NULL,mean.model=NULL,var.model=NULL,solver=0,no.col=FALSE,standard=TRUE)
oplot(out,no.col=FALSE,standard=TRUE)
JOPexampledata()
JOPexample(solver=0,optreg=0,no.col=TRUE)
locate(out,xlu=NaN,no.col=FALSE,standard=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nx}{
  number of design parameters
}
  \item{ny}{
  number of responses
}
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
  
  \code{optreg=1}: cube, 
  
  \code{optreg=0}: sphere (see Details)
}
  \item{data}{
  the data set (data.frame)
}
  \item{tau}{
  vector of target values for the corresponding responses
}
\item{interact}{
  User can decide whether the interaction effects should be included into the mean model,
  
  \code{interact=1:} interaction effects included,  
  
  \code{interaction=0:} interaction effects not included 
}
\item{quad}{
  User can decide whether the quadratic effects should be included into the mean model, 
  
  \code{quad=1:} quadratic effects included, 
  
  \code{quad=0:} quadratic effects not included 
}
\item{main.disp}{
  User can decide whether the main effects should be included into the dispersion model,
  
  \code{main.disp=1:} main effects included, 
  
  \code{main.disp=0:} main effects not included 
}
\item{interact.disp}{
  User can decide whether the interaction effects should be included into the dispersion model,
  
  \code{interact.disp=1:} interactions effects included, 
  
  \code{interact.disp=0:} interactions effects included 
}
\item{quad.disp}{
  User can decide whether the quadratic effects should be included into the dispersion model, 
  
  \code{quad.disp=1:} quadratic effects included, 
  
  \code{quad.disp=0:} quadratic effects included 
}
  \item{mean.model}{
  list of functions that model the mean for the corresponding response (see Details)
}
  \item{var.model}{
  list of functions that model the variance for the corresponding response (see Details)
}
  \item{solver}{
  optional. A numeric vector. Choices are 0, 1, 2. 
  
  If the user sets \code{solver=1},
  then JOP uses 'gosolnp' for optimization. 
  
  If the user sets \code{solver=2},
  then JOP uses 'genoud' for optimization. 
  
  Otherwise (\code{solver=0}) JOP uses nlminb for optimization. 
  
  The default is \code{solver=0}.
}
 \item{out}{object from JOP 
}
  \item{xlu}{The x-coordinate that indicates where the user assumes a "good" compromise,
  see Details  
  }
    \item{no.col}{ 
  If TRUE the plot will be gray scaled. Otherwise the plot will be coloured.
  }
    \item{standard}{ 
  If TRUE the standard deviations will be displayed on the right hand plot. 
  }
}

\value{
JOP returns a list containing the following elements:
\item{Parameters}{The i-th row of this matrix contains the optimal Parameter setting appropriate to the i-th weight matrix} 
\item{Responses}{The i-th row of this matrix contains the predicted Responses appropriate to the i-th weight matrix}  
\item{StandardDeviation}{The i-th row of this matrix contains the standard deviation value for each response}
\item{OptimalValue}{This vector contains the optimal value of the risk function for each optimal parameter setting}
\item{TargetValue}{Contains the target values for the correspoding responses}
\item{DGLM}{If no models assigned then the list DGLM contains the calculated models for the mean and dispersion for every response}
}
\details{Let \eqn{Y=(Y_1,\ldots,Y_p)}{Y=(Y_1,...,Y_p)} denote a vector of responses, \eqn{\tau=(\tau_1,\ldots,\tau_p)}{
tau=(tau_1,...,tau_p)} a vector of target values and \eqn{x=(x_1,\ldots,x_k)}{x=x_1,...,x_k} a vector of design parameters. 
 First JOP minimizes the risk function (introduced by Pignatiello (1993)) \deqn{\hat{R}(x)=\mbox{trace}(C\hat{\Sigma}(x))+
 (\hat{\mu}(x)-\tau)^TC(\hat{\mu}(x)-\tau)\mbox{,}}{R(x)=trace(C*Sigma(x)+(mu(x)-tau)'C(mu(x)-tau))} 
 where \eqn{\mbox{C}}{C} is the so called cost matrix given by \eqn{C=A^TWA}{A'WA}, where A is the so called standadization
 matrix and W is the weight matrix.
The matrix \eqn{\mbox{C}}{C} is unaffected by affine linear transformations. 
\eqn{A} is given by \deqn{A=\mbox{diag}([\frac{1}{k}\sum_{i=1}^{k}\hat{\mbox{Var}}(Y_r|x_i)]^{-\frac{1}{2}}_
{r=1,\ldots,p})\mbox{.}}{A=diag((1/k)*sum_{i=1}^{k} Var(Y_r|x_i)^(-1/2) (for r=1,...,p)).}
The weight matrices are calculated with the help of 
\eqn{\mbox{d}}{d}, \eqn{\mbox{Wstart}}{Wstart}, \eqn{\mbox{Wend}}{Wend}, \eqn{\mbox{numbW}}{numbW}. The number of weight 
matrices is equal to numbW and the weight matrices are calculated as illustrated in the following. 
Let \deqn{a=\{a_1,\ldots,a_{\mbox{numbW}}\}=\{\mbox{Wstart}, ... , \mbox{Wend}\}}{a=(a_1, ... ,a_{numbW}={Wstart, ... ,Wend}} be
a finite equidistant sequence. 
Furthermore let \eqn{\mbox{diagW}_i=(w_1,\ldots,w_p)}{diagWi=(w_1,...,w_p)}
 denote the diagonal of the i-th weight matrix \eqn{\mbox{W}_i}{Wi}. 
Then \eqn{\mbox{diagW}_i}{diagWi} is calculated in the following way: 
\deqn{\log(\mbox{diagW}_i)=d*a_i,\quad i=1\mbox{,} \ldots \mbox{,} \mbox{numbW}}{log(diagW_i)=d*a_i ,  i=1 , ... , numbW}    

The functions in the lists \code{mean.model} and \code{var.model} must depend on the parameter vector x and should only return the function value. 
If the user does not input the lists \code{mean.model} and \code{var.model} then \code{JOP} builds double generalized linear models for the mean and dispersion
for every response by means of the package \code{JointModeling}. The user has to indicate if the main, interaction or quadratic effects have to be 
included in the formula for both mean and dispersion. We assume a normal distribution with the canonical link for the mean model and 
a gamma distribution for the dispersion with log-link. Starting with the full model including all selected effects \code{JOP} conducts a backward selection 
and successively removes the least significant covariate.  

After the data is prepared, the optimization routine starts and calculates the optimal design parameters and the appropriate predicted responses
for every weight matrix and every response variable. The optimization region is a sphere around the origin if the user chooses \code{optreg=0}
and a polyhedron if the user chooses \code{optreg=1}. The radius of the sphere 
is the largest norm of the parameter vectors out of the design matrix. The polyhedron is bounded above by the largest and below by the smallest value
from each column of the matrix containing the responses. The user may choose between the optimization routines \code{nlminb} (by default)
, \code{gosolnp} (by setting \code{solver=1}, requires the package \code{Rsolnp}) and genoud (by setting \code{solver=2}, requires the package \code{rgenoud}). The routine \code{nlminb} is a local optimizer. \code{nlminb} calculates an
optimal parameter setting for three different start parameters and chooses the one with lowest objective value. Furthermore \code{nlminb} is less time-consuming then 
\code{gosolnp}. The routine \code{gosolnp} enables the generation of a set of randomly chosen parameters from which to initialize multiple
restarts of the solver \code{solnp}. \code{solnp} belongs to the class of indirect solvers and implements the augmented Lagrange multiplier
method with an SQP interior algorithm. The most time consuming solver compared to the last mentioned is \code{genoud}. This routine combines evolutionary 
search algorithms with derivative-based (Newton or quasi-Newton) methods. If the objective function is nonsmooth or has many local minima,
the user should choose the \code{gosolnp} (\code{solver=1}). Furthermore it is recommendable to choose \code{genoud} (\code{solver=2}), if the objective function
is very complex and the optimization problem is rather difficult. 

Afterwards the joint optimization plot is displayed. The left-hand plot shows the optimal design parameters for every weight matrix 
and every response variable. The legend in the left upper corner displays the names of the appropriate parameters. 
The right-hand plot shows the appropriate predicted responses. The curves correspond to the values on the vertical axis 
in the same color. The legend in the upper right corner contains the names of the appropriate responses and the target
values. On the right-hand plot the user is able to find a  "good" compromise and afterwards the user looks for the corresponding optimal 
design parameters on the left-hand plot. The dashed horizontal lines display the target values. 

The function \code{oplot} displays the joint optimization plot with the help of the object \code{out} from \code{JOP}. 

Furthermore \code{JOPexampledata} returns a list that contains the data for a sheet metal hydroforming process as an illustrating example. There 
There are two parameters and two response. The list consists of a data frame and two lists of functions that model the mean
and the variance of the corresponding response.  

The function \code{JOPexample} calculates the joint optimization plot for the example mentioned above
and returns the same object as \code{JOP}. 

The function \code{locate} allows the user to choose a preferred point on the right hand plot with a click with the mouse. 
Another possibility is to select the value \code{xlu} that has be chosen in the interval \eqn{\code{[1,\ldots,\mbox{numbW}]}}{\code{[1,...,numbW]}}. 
The value \code{xlu} indicates where
the user assumes a "good" compromise on the right hand plot. Afterwards a vertical line through the 
chosen x-coordinate is displayed on the left hand and the right hand plot. The intersection points
are the chosen responses on the right hand plot and the corresponding design parameters on the left 
hand plot. In addition the function \code{locate} returns a list that contains the chosen responses, its proper
design parameters and calculated double generalized linear models.
}
\references{
Sonja Kuhnt and Martina Erdbruegge (2004), "A strategy of robust paramater design for multiple responses", 
Statistical Modelling; 4: 249-264, TU Dortmund.

Joseph J. Pignatiello (1993), "Strategies for robust multiresponse quality engineering", IIE Transactions 25, 5-15, Texas A&M University.

Walter R. Mebane,Jr., and Jasjeet Singh Sekhon (2007), "Genetic Optimization Using Derivatives: The rgenoud package for R." Journal of Statistical Software, Forthcoming. \url{http://www.jstatsoft.org/}.  

Alexios Ghalanos and Stefan Theussl (2010). Rsolnp: General Non-linear Optimization Using Augmented Lagrange Multiplier Method. R package version 1.0-6. 

Mathieu Ribatet and Bertrand Iooss (2010). JointModeling: Joint Modelling of Mean and Dispersion. R package version 1.0-2.
}
\author{
Sonja Kuhnt and Nikolaus Rudak
}


\examples{
# Example: Sheet metal hydroforming process
data<-JOPexampledata()
# the first argument contains the data set
dataset<-data[[1]]
# JOPexample produces the same object out as JOP
# It calculates a solution for the data set mentioned above
out<-JOPexample()
# Visualization with oplot
oplot(out)
# The user may locate a point on the right plot
# and afterwards a vertical line is displayed 
# on the left and right plot.
#
# Choose the coordinate: 
xlu<-2
opt<-locate(out,xlu)
}

