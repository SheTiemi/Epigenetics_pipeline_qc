setwd("/gpfs/ysm/scratch60/montalvo-ortiz/stn7/YalePenn/QC_finalscript/Epigenetics_env/")
date = "2132023"

library(ENmix)

#Open the pheno file
pheno = read.csv("Meffil_pheno.csv")
row.names(pheno) = pheno$SampleID 
dim(pheno)

rgSet <- readidat(path="RawData",force=TRUE)
rgSet

plotCtrl(rgSet,colnames(rgSet))
#colnames indicates the order of the samples

qc<-QCinfo(rgSet)
save(rgSet,gc,pheno, fiie=paste("ENmix_qc_",date,".RData", sep=""))

#Data distribution plots
mraw <- getmeth(rgSet)
mraw

#total intensity plot is useful for data quality inspection and identification of outlier samples
pdf(file="Enmix_Meffil_distribution.pdf")
multifreqpoly(assays(mraw)$Meth+assays(mraw)$Unmeth, xlab="Total intensity")
dev.off()

beta<-getB(mraw)
anno=rowData(mraw)
beta1=beta[anno$Infinium_Design_Type=="I",]
beta2=beta[anno$Infinium_Design_Type=="II",]

library(geneplotter)
jpeg("Enmix_Meffil_dist.jpg",height=900,width=600)
par(mfrow=c(3,2))
multidensity(beta,main="Multidensity")
multifreqpoly(beta,main="Multifreqpoly",xlab="Beta value")
multidensity(beta1,main="Multidensity: Infinium I")
multifreqpoly(beta1,main="Multifreqpoly: Infinium I", xlab="Beta value")
multidensity(beta2,main="Multidensity: Infinium II")
multifreqpoly(beta2,main="Multifreqpoly: Infinium II", xlab="Beta value")
dev.off()
rm(beta1, beta2)

#filtering cross probes, polymorphic, non-specific and those detected with low quality in meffil step

info_raw = getCGinfo(rgSet)
SNPs = na.omit(c(info_raw$Probe_rs,info_raw$CpG_rs))
length(SNPs)

filter_data <- function(X){
    nonspecific_probes <- read.table("/gpfs/ysm/project/montalvo-ortiz/stn7/scripts/removalCpGs/nonspecific_probes.txt")[,1]
    polymorphic_probes <- read.table("/gpfs/ysm/project/montalvo-ortiz/stn7/scripts/removalCpGs/polymorphic_cpgs.txt")[,1]
    cross = read.csv("/gpfs/ysm/project/montalvo-ortiz/stn7/scripts/removalCpGs/cross.csv",stringsAsFactors=FALSE,header=TRUE)[,1]
    meffil = read.csv("CpGremoval/Meffil_badCpGs.csv",stringsAsFactors=FALSE,header=TRUE)
    exclude <- union(union(nonspecific_probes,union(polymorphic_probes,cross)),meffil)
    return(X[setdiff(rownames(X),exclude),])
}

`%notin%` <- Negate(`%in%`)
qc.beta <- filter_data(beta)
exclude = row.names(mraw[row.names(mraw) %notin% row.names(qc.beta)])
length(exclude)

exclude = union(exclude,SNPs)
length(exclude)

#background correction and dye bias correction
#if qc info object is provided, the low quality or outlier samples and probes will be excluded before background correction
mdat<-preprocessENmix(rgSet, bgParaEst="oob", dyeCorr="RELIC", QCinfo=qc, exQCsample=TRUE, exQCcpg=TRUE, nCores=6, exCpG=exclude)

#between-array normalization
#"quantile1", "quantile2", or "quantile3".
mdat_q1<-norm.quantile(mdat, method="quantile1")

#probe-type bias adjustment
beta_q1<-rcp(mdat_q1,qcscore=qc)

#batch effect correction
sva<-ctrlsva(rgSet)

save(mdat,mdat_q1,beta_q1,sva, file=paste("ENmix_normalized_quantile1_",date,".RData", sep=""))

system("rm -r RawData")
