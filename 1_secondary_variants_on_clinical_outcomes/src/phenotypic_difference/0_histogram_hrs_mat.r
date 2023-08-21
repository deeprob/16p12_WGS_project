#Note: HRS-MAT scores are for probands/carriers/noncarriers in integration table
mydata<-read.table("hrs_mat_scores_0822.txt",sep="\t",header=T)

hist(mydata$proband, breaks=seq(0,200,l=100), xlim=c(0,200), ylim=c(0,0.125), col=rgb(0,0,1,0.5), xlab="HRS-MAT",freq=F)
hist(mydata$carrier_parent, breaks=seq(0,200,l=100), col=rgb(1,0,0,0.5),add=T,freq=F)
hist(mydata$noncarrier_parent, breaks=seq(0,200,l=100), col=rgb(0,0.545,0,0.5), add=T,freq=F)

lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$proband)), sd=sqrt(var(na.omit(mydata$proband)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$carrier_parent)), sd=sqrt(var(na.omit(mydata$carrier_parent)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$noncarrier_parent)), sd=sqrt(var(na.omit(mydata$noncarrier_parent)))), col="green4",lwd=4)

#Plot histogram lines only
plot(1, type="n", xlab="HRS-MAT", ylab="", xlim=c(0, 200), ylim=c(0, 0.04))
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$proband)), sd=sqrt(var(na.omit(mydata$proband)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$carrier_parent)), sd=sqrt(var(na.omit(mydata$carrier_parent)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$noncarrier_parent)), sd=sqrt(var(na.omit(mydata$noncarrier_parent)))), col="green4",lwd=4)

#Histogram comparing 16p12.1 and SSC cohort (SSC values from Hanson, JADD 2014)
plot(1, type="n", xlab="HRS-MAT score", ylab="", xlim=c(0, 200), ylim=c(0, 0.04))
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=83, sd=15), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=99, sd=18), col="red",lwd=4)
