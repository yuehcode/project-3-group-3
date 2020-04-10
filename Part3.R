

# load the 9 treated sampls.
data1 <- read.table("SRR1177981Aligned.sortedByCoord.out.txt", header=TRUE, skip=1)
data2<- read.table("SRR1177982Aligned.sortedByCoord.out.txt",header=TRUE, skip=1)
data3<- read.table("SRR1177983Aligned.sortedByCoord.out.txt",header = TRUE,skip=1)
data4<- read.table("SRR1178008Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)
data5<- read.table("SRR1178009Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)
data6<-read.table("SRR1178010Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)
data7<-read.table("SRR1178014Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)
data8<-read.table("SRR1178021Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)
data9<-read.table("SRR1178047Aligned.sortedByCoord.out.txt",header=TRUE,skip=1)

#bind 9 treated sampls to one data file
data<-cbind(data1[7],data2[7],data3[7],data4[7],data5[7],data6[7],data7[7],data8[7],data9[7])

#make column name
sampleNames <- c("SRR1177981", "SRR1177982", "SRR1177983", "SRR1178008", 
"SRR1178009", "SRR1178010","SRR1178014","SRR1178021","SRR1178047")
names(data)[1:9] <- sampleNames

#add the gene name
data<- cbind(data1[1],data)

##save to one file
write.csv(data,file="concatenated_counts.csv")

# Get log2 counts per million and make boxplot to show distribution!
# Check distributions of samples using boxplots
par(mfrow=c(1,1))
logcounts <- log2(data[2:10] + 1)
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,cex.axis=0.8)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")

# load control file
controldata<- read.csv("control_counts.csv", header = T)
# make the name for the controls 
corn_oil <- c("SRR1178050","SRR1178061","SRR1178063")
saline <- c("SRR1178004","SRR1178006","SRR1178013")

# make the name list for treated samples into each group!
AhR_corn_oil <- c("SRR1178008","SRR1178009", "SRR1178010")
CAR_corn_oil <- c("SRR1178014","SRR1178021","SRR1178047")

dna_damage_saline <- c("SRR1177981", "SRR1177982", "SRR1177983")

# subset to comparison file (dna_damaged treated vs samples)
dna_damage_saline_deseq<- cbind(data1[1],data[,dna_damage_saline], controldata[, saline ])
# save the count file for later use
write.csv(dna_damage_saline_deseq,file="dna_damage_saline_deseq.csv",row.names = F)

# subset to  each comparison file (dna_damaged treated vs samples)
AhR_corn_oil_deseq <- cbind(data1[1],data[,AhR_corn_oil], controldata[, corn_oil])
# save the count file for later use
write.csv(AhR_corn_oil_deseq,file="AhR_corn_oil_deseq.csv", row.names = F)


#  subset to comparison file (dna_damaged treated vs samples)
CAR_corn_oil_deseq <- cbind(data1[1],data[,CAR_corn_oil], controldata[, corn_oil])
# save the count file for later use
write.csv(CAR_corn_oil_deseq,file="CAR_corn_oil_deseq.csv", row.names = F)



