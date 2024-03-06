n<-32
# we want to keep the mean remain the same 
#which is zeta_shape x zeta_scale is remain the same
zeta_shape<-128   
zeta_scale<-10/zeta_shape
zeta_m<- rgamma(32, shape=zeta_shape,scale= zeta_scale)+1

A<-seq(0,50,by=0.1)
#dgamma is the density function 
#not sure why we have A-1
#forgot it
plot(x=A,y=dgamma(A-1,shape = zeta_shape, rate=1/zeta_scale),type="l")
#another question in the actual code have to make sure it is the one we choosed 
#here just random testing we didn't use it in the simulation yet
#so we either need to plot it each time? or we just leave it as test just plot the entropy 

points(x=zeta_m,y=rep(0,n),xlim=c(0,50),col="red")
