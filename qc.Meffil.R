setwd("/gpfs/ysm/scratch60/montalvo-ortiz/stn7/YalePenn/QC_finalscript/Epigenetics_env/")
dataDirectory <- (".")
list.files(dataDirectory, recursive = TRUE)

date="2132023"

########################  ###############################

library(meffil)

dir = read.csv("EWastools_pheno.csv")
row.names(dir) = dir$SampleID
dim(dir)

colnames(dir)[2] <- "Basename"
dir$Sex = dir$A1_SEX
dir$Sex[dir$A1_SEX==1] = "M"
dir$Sex[dir$A1_SEX==2] = "F"

#Create samplesheet
samplesheet = meffil.create.samplesheet(dir$Basename,basename=dir$Basename)
samplesheet$Sex = dir$Sex[match(dir$SampleID,samplesheet$Sample_Name)]

#Run the qc
qc.objects = meffil.qc(samplesheet, cell.type.reference = "blood gse35069 complete", verbose = TRUE)
#Create summary
qc.summary = meffil.qc.summary(qc.objects, verbose = TRUE)
#saving the badCpGs
system("mkdir CpGremoval")
write.csv(qc.summary$bad.cpgs, file="CpGremoval/Meffil_badCpGs.csv",row.names=FALSE)

#Create report
meffil.qc.report(qc.summary,
                 output.file = paste("Meffil_qc_",date,".html",sep=""),
                 author = "Sheila",
                 study = "QC") 

#Meffil plots after removing outliers
paste("remove outlier samples if necessary")
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

pdf(file=paste("Meffilplotpc_",date,".pdf",sep=""))
print(meffil.plot.pc.fit(qc.objects)$plot)
dev.off()

#Change the dir file
dir$Meffil_QC = FALSE
dir$Meffil_QC[dir$SampleID %in% qc.summary$bad.samples$sample.name] = TRUE
sum(dir$Meffil_QC)
dir = dir[dir$Meffil_QC==FALSE,]

write.csv(dir, "Meffil_pheno.csv",row.names=FALSE)

#Creating the RawData directory with the link for all files to be analyzed with Enmix
x = dir$Basename

system("mkdir RawData")
write.csv(x, file="RawData/directions.csv",row.names=FALSE)
cdm = paste("perl /gpfs/ysm/project/montalvo-ortiz/stn7/scripts/ln.pl","RawData/directions.csv","RawData/.",sep=" ")
system(cdm)

