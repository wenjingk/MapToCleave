dumbbell_heat_plot<-function(ids,title){

#ids<-ids_inducible
#title<-"test"

# load bulge profile
all_ids<-ids
dir<-"./MapToCleave_supplemental_data/detailed_feature_profile_of_MapToCleave_processed_and_unprocessed_miRNA_precursors/"

bulge_5p<-vector()

for (i in 1:length(all_ids)){
    file<-paste(paste(dir,all_ids[i],sep="/"),".HEK.unopened_5p_ref",sep="")
    # If file is not exist, this indicates the structure is not folded as long hairpin. Assum open structure *****
    if (!file.exists(file)){
        #tmp<-data.frame(seq(-25,30),rep("*",times = length(seq(-25,30))))
        next
    } else {
        tmp<-read.table(file,sep="\t")
    }
    colnames(tmp)<-c("position",all_ids[i])
    if (i == 1){
        bulge_5p<-tmp
    } else {
        bulge_5p<-merge(bulge_5p,tmp,x.by="position",y.by="position",all=TRUE)
    }
}

#bulge_5p<-bulge_5p[order(bulge_5p$position),]

for (i in 2:length(bulge_5p[1,])){
    bulge_5p[,i]<-as.character(bulge_5p[,i])
}

bulge_5p<-bulge_5p[bulge_5p$position %in% seq(-20,30),]
rownames(bulge_5p)<-bulge_5p[,1]
bulge_5p[,1]<-NULL


bulge_3p<-vector()

for (i in 1:length(all_ids)){
    file<-paste(paste(dir,all_ids[i],sep="/"),".HEK.unopened_3p_ref",sep="")
    # If file is not exist, this indicates the structure is not folded as long hairpin. Assum open structure *****
    if (!file.exists(file)){
        #tmp<-data.frame(seq(-25,30),rep("*",times = length(seq(-25,30))))
        next
    } else {
        tmp<-read.table(file,sep="\t")
    }
    colnames(tmp)<-c("position",all_ids[i])
    if (i == 1){
        bulge_3p<-tmp
    } else {
        bulge_3p<-merge(bulge_3p,tmp,x.by="position",y.by="position",all=TRUE)
    }
}

#bulge_3p<-bulge_3p[order(bulge_3p$position),]

for (i in 2:length(bulge_3p[1,])){
    bulge_3p[,i]<-as.character(bulge_3p[,i])
}

bulge_3p<-bulge_3p[bulge_3p$position %in% seq(-20,30),]
rownames(bulge_3p)<-bulge_3p[,1]
bulge_3p[,1]<-NULL


# load bulge features
all_ids<-ids
dir<-"./MapToCleave_supplemental_data/detailed_feature_profile_of_MapToCleave_processed_and_unprocessed_miRNA_precursors/"

bulge_features_5p<-vector()

for (i in 1:length(all_ids)){
    file<-paste(paste(dir,all_ids[i],sep="/"),".HEK.bulgefeatures_unopened_5p_ref",sep="")
    if (!file.exists(file)){
        next
        } else {
        tmp<-read.table(file,sep="\t")
        bulge_features_5p<-rbind(bulge_features_5p,tmp)
    }
}

colnames(bulge_features_5p)<-c("ID","not-important","bulge_size","bulge_start","bulge_end","bulge_type","bulge_5pseq","bulge_3pseq","bulge_5pend","bulge_3pend")
bulge_features_5p[,1]<-as.character(bulge_features_5p[,1])
bulge_features_5p[,1]<-gsub(" mirtron","",bulge_features_5p[,1])


bulge_features_3p<-vector()

for (i in 1:length(all_ids)){
    file<-paste(paste(dir,all_ids[i],sep="/"),".HEK.bulgefeatures_unopened_3p_ref",sep="")
    if (!file.exists(file)){
        next
        } else {
        tmp<-read.table(file,sep="\t")
        bulge_features_3p<-rbind(bulge_features_3p,tmp)
    }
}

colnames(bulge_features_3p)<-c("ID","not-important","bulge_size","bulge_start","bulge_end","bulge_type","bulge_5pseq","bulge_3pseq","bulge_5pend","bulge_3pend")
bulge_features_3p[,1]<-as.character(bulge_features_3p[,1])
bulge_features_3p[,1]<-gsub(" mirtron","",bulge_features_3p[,1])

## bulge_size_plot 

# get information when using 5p unopened stem as coordinates
bulge_features_subset<-bulge_features_5p[bulge_features_5p$ID %in% ids,]
features_subset<-t(bulge_5p[,colnames(bulge_5p) %in% ids,drop=F])
features_subset<-as.data.frame(features_subset)
    
     
# get bulge percentage
bulge_percent<-vector()
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == ">" | features_subset[,j] == "<" | features_subset[,j] == "*"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(as.numeric(colnames(features_subset)),tmp)
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == "*"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(bulge_percent,tmp)
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == ">" | features_subset[,j] == "<"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(bulge_percent,tmp)


# get bulge size distribution
bulge_features_percentage<-vector()
tmp<-vector()
nonNA_and_nonparing_count<-vector()
for (j in seq(-20,30)){
    tmp<-bulge_features_subset[bulge_features_subset$bulge_start <= j & bulge_features_subset$bulge_end >= j,]
    tmp<-data.frame(table(tmp$bulge_size))
    if (length(tmp[,1]) == 0){
        tmp<-data.frame(c(1,2),c(0,0))
        tmp[,1]<-as.factor(tmp[,1])
    }
    colnames(tmp)<-c("bulge_size",j)
    if (j == -20){
        bulge_features_percentage <- tmp
    } else {
        bulge_features_percentage <- merge(bulge_features_percentage,tmp,x.by="bulge_size",y.by="bulge_size",all=TRUE)
    }
    nonNA_and_nonparing_count<-c(nonNA_and_nonparing_count,sum(tmp[,2]))
}

bulge_features_percentage[is.na(bulge_features_percentage)]<-0
bulge_features_percentage[,1]<-as.numeric(as.character(bulge_features_percentage[,1]))
bulge_features_percentage<-bulge_features_percentage[order(bulge_features_percentage[,1],decreasing = FALSE),]  
bulge_ge5<-apply(bulge_features_percentage[bulge_features_percentage[,1] >= 5,2:length(bulge_features_percentage[1,])],2,sum)
bulge_features_percentage<-bulge_features_percentage[bulge_features_percentage[,1] < 5,]
bulge_features_percentage[,1]<-as.character(bulge_features_percentage[,1])
bulge_features_percentage<-rbind(bulge_features_percentage,c("ge5",bulge_ge5))
rownames(bulge_features_percentage)<-bulge_features_percentage[,1]
bulge_features_percentage[,1]<-NULL

# nonNA_count<-vector()
# for (j in 1:length(features_subset[1,])){
#     nonNA_count<-c(nonNA_count,length(which(!is.na(features_subset[,j]))))
# }

# for ( i in 1:length(bulge_features_percentage[1,])){
#     bulge_features_percentage[,i]<-as.numeric(as.character(bulge_features_percentage[,i]))
#     bulge_features_percentage[,i] <- bulge_features_percentage[,i]/nonNA_count[i]*100
# }

for ( i in 1:length(bulge_features_percentage[1,])){
    bulge_features_percentage[,i]<-as.numeric(as.character(bulge_features_percentage[,i]))
    bulge_features_percentage[,i] <- bulge_features_percentage[,i]/nonNA_and_nonparing_count[i]*bulge_percent[i,2]
}

bulge_features_percentage_5p<-bulge_features_percentage
bulge_percent_5p<-bulge_percent
bulge_features_subset_5p<-features_subset

# get information when using 3p unopened stem as coordinates
bulge_features_subset<-bulge_features_3p[bulge_features_3p$ID %in% ids,]
features_subset<-t(bulge_3p[,colnames(bulge_3p) %in% ids,drop=F])
features_subset<-as.data.frame(features_subset)
    
     
# get bulge percentage
bulge_percent<-vector()
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == ">" | features_subset[,j] == "<" | features_subset[,j] == "*"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(as.numeric(colnames(features_subset)),tmp)
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == "*"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(bulge_percent,tmp)
tmp<-vector()
for (j in 1:length(features_subset[1,])){
    tmp<-c(tmp,length(which(features_subset[,j] == ">" | features_subset[,j] == "<"))/length(which(!is.na(features_subset[,j])))*100)
}
bulge_percent<-cbind(bulge_percent,tmp)

# get bulge size distribution
bulge_features_percentage<-vector()
tmp<-vector()
nonNA_and_nonparing_count<-vector()
for (j in seq(-20,30)){
    tmp<-bulge_features_subset[bulge_features_subset$bulge_start <= j & bulge_features_subset$bulge_end >= j,]
    tmp<-data.frame(table(tmp$bulge_size))
    if (length(tmp[,1]) == 0){
        tmp<-data.frame(c(1,2),c(0,0))
        tmp[,1]<-as.factor(tmp[,1])
    }
    colnames(tmp)<-c("bulge_size",j)
    if (j == -20){
        bulge_features_percentage <- tmp
    } else {
        bulge_features_percentage <- merge(bulge_features_percentage,tmp,x.by="bulge_size",y.by="bulge_size",all=TRUE)
    }
    nonNA_and_nonparing_count<-c(nonNA_and_nonparing_count,sum(tmp[,2]))
    
}

bulge_features_percentage[is.na(bulge_features_percentage)]<-0
bulge_features_percentage[,1]<-as.numeric(as.character(bulge_features_percentage[,1]))
bulge_features_percentage<-bulge_features_percentage[order(bulge_features_percentage[,1],decreasing = FALSE),]  
bulge_ge5<-apply(bulge_features_percentage[bulge_features_percentage[,1] >= 5,2:length(bulge_features_percentage[1,])],2,sum)
bulge_features_percentage<-bulge_features_percentage[bulge_features_percentage[,1] < 5,]
bulge_features_percentage[,1]<-as.character(bulge_features_percentage[,1])
bulge_features_percentage<-rbind(bulge_features_percentage,c("ge5",bulge_ge5))
rownames(bulge_features_percentage)<-bulge_features_percentage[,1]
bulge_features_percentage[,1]<-NULL

# nonNA_count<-vector()
# for (j in 1:length(features_subset[1,])){
#     nonNA_count<-c(nonNA_count,length(which(!is.na(features_subset[,j]))))
# }

# for ( i in 1:length(bulge_features_percentage[1,])){
#     bulge_features_percentage[,i]<-as.numeric(as.character(bulge_features_percentage[,i]))
#     bulge_features_percentage[,i] <- bulge_features_percentage[,i]/nonNA_count[i]*100
# }


for ( i in 1:length(bulge_features_percentage[1,])){
    bulge_features_percentage[,i]<-as.numeric(as.character(bulge_features_percentage[,i]))
    bulge_features_percentage[,i] <- bulge_features_percentage[,i]/nonNA_and_nonparing_count[i]*bulge_percent[i,2]
}


bulge_features_percentage_3p<-bulge_features_percentage
bulge_percent_3p<-bulge_percent
bulge_features_subset_3p<-features_subset

# get pair percentage
features_subset_5p<-bulge_features_subset_5p
features_subset_3p<-bulge_features_subset_3p

pair_percent_5p<-vector()
tmp<-vector()
for (j in 1:length(features_subset_5p[1,])){
  tmp<-c(tmp,length(features_subset_5p[which(features_subset_5p[,j] == "C-G" | features_subset_5p[,j] == "G-C"),,drop=F][,1])/length(which(!is.na(features_subset_5p[,j])))*100)
}
pair_percent_5p<-cbind(as.numeric(colnames(features_subset_5p)),tmp)
tmp<-vector()
for (j in 1:length(features_subset_5p[1,])){
  tmp<-c(tmp,length(features_subset_5p[which(features_subset_5p[,j] == "A-T" | features_subset_5p[,j] == "T-A"),,drop=F][,1])/length(which(!is.na(features_subset_5p[,j])))*100)
}
pair_percent_5p<-cbind(pair_percent_5p,tmp)
tmp<-vector()
for (j in 1:length(features_subset_5p[1,])){
  tmp<-c(tmp,length(features_subset_5p[which(features_subset_5p[,j] == "G-T" | features_subset_5p[,j] == "T-G"),,drop=F][,1])/length(which(!is.na(features_subset_5p[,j])))*100)
}
pair_percent_5p<-cbind(pair_percent_5p,tmp)
pair_percent_5p<-data.frame(pair_percent_5p)
rownames(pair_percent_5p)<-pair_percent_5p[,1]
pair_percent_5p[,1]<-NULL
colnames(pair_percent_5p)<-c("GC","AU","GU")
pair_percent_5p<-pair_percent_5p[,c("GU","AU","GC")]
pair_percent_5p<-t(pair_percent_5p)

pair_percent_3p<-vector()
tmp<-vector()
for (j in 1:length(features_subset_3p[1,])){
  tmp<-c(tmp,length(features_subset_3p[which(features_subset_3p[,j] == "C-G" | features_subset_3p[,j] == "G-C"),,drop=F][,1])/length(which(!is.na(features_subset_3p[,j])))*100)
}
pair_percent_3p<-cbind(as.numeric(colnames(features_subset_3p)),tmp)
tmp<-vector()
for (j in 1:length(features_subset_3p[1,])){
  tmp<-c(tmp,length(features_subset_3p[which(features_subset_3p[,j] == "A-T" | features_subset_3p[,j] == "T-A"),,drop=F][,1])/length(which(!is.na(features_subset_3p[,j])))*100)
}
pair_percent_3p<-cbind(pair_percent_3p,tmp)
tmp<-vector()
for (j in 1:length(features_subset_3p[1,])){
  tmp<-c(tmp,length(features_subset_3p[which(features_subset_3p[,j] == "G-T" | features_subset_3p[,j] == "T-G"),,drop=F][,1])/length(which(!is.na(features_subset_3p[,j])))*100)
}
pair_percent_3p<-cbind(pair_percent_3p,tmp)
pair_percent_3p<-data.frame(pair_percent_3p)
rownames(pair_percent_3p)<-pair_percent_3p[,1]
pair_percent_3p[,1]<-NULL
colnames(pair_percent_3p)<-c("GC","AU","GU")
pair_percent_3p<-pair_percent_3p[,c("GU","AU","GC")]
pair_percent_3p<-t(pair_percent_3p)

number<-paste("n = ",length(unique(rownames(features_subset_5p))),sep="")

bulge_pair_df_5p<-rbind(as.matrix(bulge_features_percentage_5p),pair_percent_5p)
bulge_pair_df_5p[is.na(bulge_pair_df_5p)]<-0
bulge_pair_df_3p<-rbind(as.matrix(bulge_features_percentage_3p),pair_percent_3p)
bulge_pair_df_3p[is.na(bulge_pair_df_3p)]<-0



lty.o <- par("lty")
par(lty = 0)
cols<-c(rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090')),rev(c("#2166ac","#92c5de","#f7f7f7")))
col_names<-c("1","2","3","4","ge5","GU","AU","GC")
col_df<-data.frame(col_names,cols)
col_df[,1]<-as.character(col_df[,1])
col_df[,2]<-as.character(col_df[,2])
cols_5p<-col_df[col_df[,1] %in% rownames(bulge_pair_df_5p),][,2]


df.bar<-barplot(bulge_pair_df_5p,
                col=cols_5p,
                legend.text = FALSE,ylim=c(-105,105),las=2,space=0, cex.names=0.8,
                ylab="Bulge [%] & bulge size", xlab="nt from Drosha cutting site",
                main=title,xaxt='n')

cols_3p<-col_df[col_df[,1] %in% rownames(bulge_pair_df_3p),][,2]


df.bar<-barplot(bulge_pair_df_3p*-1,
                col=cols_3p,
                legend.text = FALSE,ylim=c(-105,105),las=2,space=0, cex.names=0.8,add = TRUE,xaxt='n')

par(lty = lty.o)
points(df.bar,bulge_percent_5p[,2],pch=19,cex=0.5)
lines(df.bar,bulge_percent_5p[,2],lwd = 0.8)
points(df.bar,bulge_percent_3p[,2]*-1,pch=19,cex=0.5)
lines(df.bar,bulge_percent_3p[,2]*-1,lwd = 0.8)
abline(h=0,lty=1)
abline(v=7,lty=2,col="white")
abline(v=14,lty=2,col="white")
abline(v=20,lty=1,col="black")
abline(v=42,lty=2,col="white")
x_axis<-colnames(as.matrix(bulge_features_percentage_5p))
axis(1, at = df.bar[seq(1, length(df.bar), 2)], las=2, labels = x_axis[seq(1, length(x_axis), 2)],tick = TRUE)
legend("top",legend=number,bty = "n")

}
