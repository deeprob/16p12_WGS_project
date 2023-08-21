#Note: SRS scores are for probands/carriers/noncarriers in integration table
mydata<-read.table("srs_scores_1022.txt",sep="\t",header=T)

hist(mydata$proband, breaks=seq(0,200,l=100), xlim=c(0,200), ylim=c(0,0.06), col=rgb(0,0,1,0.5), xlab="SRS",freq=F)
hist(mydata$carrier_parent, breaks=seq(0,200,l=100), col=rgb(1,0,0,0.5),add=T,freq=F)
hist(mydata$noncarrier_parent, breaks=seq(0,200,l=100), col=rgb(0,0.545,0,0.5), add=T,freq=F)

lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$proband)), sd=sqrt(var(na.omit(mydata$proband)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$carrier_parent)), sd=sqrt(var(na.omit(mydata$carrier_parent)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$noncarrier_parent)), sd=sqrt(var(na.omit(mydata$noncarrier_parent)))), col="green4",lwd=4)

#Plot histogram lines only
plot(1, type="n", xlab="SRS", ylab="", xlim=c(0,200), ylim=c(0, 0.02))
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$proband)), sd=sqrt(var(na.omit(mydata$proband)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$carrier_parent)), sd=sqrt(var(na.omit(mydata$carrier_parent)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$noncarrier_parent)), sd=sqrt(var(na.omit(mydata$noncarrier_parent)))), col="green4",lwd=4)

#Comparison between SSC, SVIP, and our cohort
#Note this uses raw scores, not t-scores
mydata<-read.table("srs_16p12_SSC_SVIP.txt",sep="\t",header=T)

hist(mydata$X16p12.1_deletion, breaks=seq(0,200,l=100), xlim=c(0,200), ylim=c(0,0.04), col=rgb(0,0,1,0.5), xlab="SRS",freq=F)
hist(mydata$SSC_cohort, breaks=seq(0,200,l=100), col=rgb(1,0,0,0.5),add=T,freq=F)
hist(mydata$SVIP_16p11.2_cohort, breaks=seq(0,200,l=100), col=rgb(0,0.545,0,0.5), add=T,freq=F)

lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$X16p12.1_deletion)), sd=sqrt(var(na.omit(mydata$X16p12.1_deletion)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$SSC_cohort)), sd=sqrt(var(na.omit(mydata$SSC_cohort)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$SVIP_16p11.2_cohort)), sd=sqrt(var(na.omit(mydata$SVIP_16p11.2_cohort)))), col="green4",lwd=4)

#Plot histogram lines only
plot(1, type="n", xlab="SRS", ylab="", xlim=c(0, 200), ylim=c(0, 0.02))
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$X16p12.1_deletion)), sd=sqrt(var(na.omit(mydata$X16p12.1_deletion)))), col="blue",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$SSC_cohort)), sd=sqrt(var(na.omit(mydata$SSC_cohort)))), col="red",lwd=4)
lines(seq(0,200,by=0.01),dnorm(seq(0,200,by=0.01), mean=mean(na.omit(mydata$SVIP_16p11.2_cohort)), sd=sqrt(var(na.omit(mydata$SVIP_16p11.2_cohort)))), col="green4",lwd=4)