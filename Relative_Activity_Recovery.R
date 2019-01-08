#####RELATIVE ACTIVITY RECOVERY###########
##### Author - Dr.Jenni Inkinen###########
#File 1needed: Relative abundance (%) OTU table with DNA and RNA samples#
#File2 needed: Pairs list file where one sample row contain column for DNA and RNA#
#start
DNAvector <- Pairs_list$DNA
RNAvector <- Pairs_list$RNA
samplecode <- Pairs_list$Code2
otutable <- otutable_abundance
indeksitDNA <- match(DNAvector, colnames((otutable)))   
indeksitRNA <- match(RNAvector, colnames((otutable)))

#create new column to "Pairs_list" to match DNA/RNA
missing.data <- is.na(indeksitRNA) + is.na(indeksitDNA)
Pairs_list[["missing"]] <- missing.data > 0
DNAvector2 <- Pairs_list[!Pairs_list[["missing"]], "DNA"]
RNAvector2 <- Pairs_list[!Pairs_list[["missing"]], "RNA"]
indeksitDNA2 <- match(DNAvector2$DNA, colnames((otutable)))   
indeksitRNA2 <- match(RNAvector2$RNA, colnames((otutable)))


#create new file and rename 
activity.dataframe <- otutable[,indeksitRNA2]-otutable[,indeksitDNA2]
colnames(activity.dataframe) <- samplecode[!Pairs_list[["missing"]]]


#print missing values
samplecode[Pairs_list[["missing"]]]


#add OTUs and taxonomy information
activity.dataframe$OTUID <- otutable$OTUID
activity.dataframe$taxonomy <- otutable$taxonomy


#save file
write.csv2(activity.dataframe, file="activity.otutable.csv")
