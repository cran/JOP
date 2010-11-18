JOPexample<-function(solver=0,optreg=0)
{
# Example: Sheet metal hydroforming process
dataset<-as.data.frame(t(matrix(c(1.00,1.00,7.9966,0.067,1.00,1.00,20.3601,0.083,
1.00,1.00,19.1287,0.067,1.00,1.00,31.0451,0.081,1.00,-1.00,38.6362,0.061,
1.00,-1.00,38.6033,0.057,1.00,-1.00,36.7798,0.039,1.00,-1.00,44.5056,0.038,
-1.00,1.00,0.9748,0.061,-1.00,1.00,9.1252,0.078,-1.00,1.00,8.5072,0.054,
-1.00,1.00,21.3566,0.063,-1.00,-1.00,33.6568,0.063,-1.00,-1.00,37.6792,0.064,
-1.00,-1.00,35.1122,0.047,-1.00,-1.00,43.9494,0.043,
-1.41,0.00,7.2792,0.054,
-1.41,0.00,21.7991,0.074,-1.41,0.00,19.1302,0.050,-1.41,0.00,36.3412,0.063,
1.41,0.00,23.8342,0.066,1.41,0.00,31.6159,0.080,1.41,0.00,35.0138,0.065,
1.41,0.00,36.5504,0.055,
0.00,1.41,1.3413,0.065,0.00,1.41,10.7754,0.092,0.00,1.41,9.0730,0.058,0.00,1.41,21.8355,0.080,
0.00,0.00,15.4458,0.066,0.00,0.00,30.3393,0.070,0.00,0.00,29.7673,0.057,0.00,0.00,35.1129,0.062,
0.00,-1.41,35.0273,0.052,0.00,-1.41,45.1586,0.046,0.00,-1.41,46.9272,0.030,
0.00,-1.41,42.3516,0.033), nrow=4)))
dimnames(dataset)[[2]] <- c("K", "D","Area","RBT")

# Number of Design Parameters:
nx<-2

# Number of Response Variables:
ny<-2

# Target Values:
tau<-c(0,0.05)


# Set the Values for the Weight Matrices
Wstart<--9.21
numbW<-11
d<-c(1,0)
Wend<-9.21

# JOP calculates the optimal design parameters and the appropriate predicted responses.
# Furthermore it produces the joint optimization plot
out<-JOP(nx=nx,ny=ny,Wstart=Wstart,Wend=Wend,numbW=numbW,d=d,optreg=optreg,tau=tau,interact=1,
                quad=1,main.disp=1,interact.disp=1,quad.disp=1,data=dataset,solver=solver)

return(out)
}