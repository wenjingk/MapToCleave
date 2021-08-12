# the script is same as the Rice shannon entropy paper

args<-commandArgs(TRUE)

dir<-"./output/"
a<-read.table(paste(dir,args[1],"_base_pair_probability_from_RNAstructure.txt",sep=""),sep="\t",skip=1,header=TRUE)
a$logn<-a[,3]*10^(-1*a[,3])

shannon<-vector()
for (i in unique(c(a[,1],a[,2]))){
  tmp<-a[a[,1] == i | a[,2] == i,]
  shannon<-rbind(shannon,c(i,sum(tmp[,4])))
}

shannon<-data.frame(shannon)
colnames(shannon)<-c("position","positional_entropy")

dir2<-"./output/"
tmp<-read.table(paste(dir2,args[1],".HEK.statistics",sep=""),sep="\t")
start_5p<-as.numeric(as.character(tmp[tmp[,1]=="start_5p",][,2])) + 1
shannon$strand_5p<-shannon[,1]-start_5p

start_3p<-as.numeric(as.character(tmp[tmp[,1]=="start_3p",][,2])) + 1
shannon$strand_3p<-shannon[,1]-start_3p+100

out_strand_5p<-shannon[shannon$strand_5p %in% seq(-20,30),c("strand_5p","positional_entropy")]
colnames(out_strand_5p)<-c("position","positional_entropy")

out_strand_3p<-shannon[shannon$strand_3p %in% seq(90,145),c("strand_3p","positional_entropy")]
colnames(out_strand_3p)<-c("position","positional_entropy")

out<-rbind(out_strand_5p,out_strand_3p)

#plot(out[,1],out[,2],pch="",ylim=c(0,2))
#lines(out[,1],out[,2])

write.table(out,paste(dir,args[1],"_positional_entropy_drosha_site_as_rice_paper.txt",sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote = FALSE)


