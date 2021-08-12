local_stability_plot<-function(input,top,tail,title){

library(gplots)
library(RColorBrewer)
library(ggplot2)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

df<-input

top_in<-top[top %in% colnames(df)]
test1<-df[,c(top_in)]
test1[is.na(test1)]<-0

tail_in<-tail[tail %in% colnames(df)]
test2<-df[,c(tail_in)]
test2[is.na(test2)]<-0

test1_mean<-apply(test1,1,mean)
test2_mean<-apply(test2,1,mean)
y_max<-max(c(test1_mean,test2_mean))
y_min<-min(c(test1_mean,test2_mean))
    
plot(names(test1_mean),test1_mean,pch=19,col="red",ylab="mean free energy", xlab="position from Drosha site",
    main=title,ylim=c(y_min-1,y_max+1))
lines(names(test1_mean),test1_mean,col="red")
points(names(test2_mean),test2_mean,col="blue",pch=19)
lines(names(test2_mean),test2_mean,col="blue")
abline(v=-13,lty=2)
legend("bottom",legend=c(paste(deparse(substitute(top))," (n = ",length(top_in)," )",sep=""),paste(deparse(substitute(tail))," (n = ",length(tail_in)," )",sep="")),col=c("red","blue"),lty=c(1,1))

}

