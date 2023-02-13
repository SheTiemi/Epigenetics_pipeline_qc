setwd("/gpfs/ysm/scratch60/montalvo-ortiz/stn7/YalePenn/QC_finalscript/Epigenetics_env/")
dataDirectory <- (".")
list.files(dataDirectory, recursive = TRUE)

date="292023"

######################## EWAStools ###############################

library(stringi)
library(magrittr)
library(data.table)
library(svd)
library(ewastools)
library(purrr)

pheno = read.csv("Manifest_dir_1302023_YalePenn_initialPhen.csv")
row.names(pheno) = pheno$SampleID

#create the meth file
meth = read_idats(pheno$Dir)

pheno = as.data.table(pheno)
pheno[,exclude:=FALSE]
meth %>% control_metrics %>% sample_failure -> pheno$failed
fails <- subset(pheno, failed == "TRUE")
message("Failed Samples: ", nrow(fails))

meth = ewastools::detectionP(meth)

# Check Dye-bias
meth %>% dont_normalize                      -> with_bias
meth %>% correct_dye_bias %>% dont_normalize -> corrected

snps = meth$manifest[probe_type=="rs" & channel=="Both"]$index

pdf(file=paste("Ewastools_density_bias_",date,".pdf",sep=""))
plot (density(with_bias[snps,14],na.rm=TRUE,bw=0.1),col=1)
lines(density(corrected[snps,14],na.rm=TRUE,bw=0.1),col=2)
abline(v=0.5,lty=3)
legend("topleft",col=1:2,legend=c("with_bias","corrected"),lwd=1)
dev.off()

meth %>% correct_dye_bias -> meth
rm(corrected,with_bias)
gc()

meth = mask(meth,0.01)

#Sex mismatches
chrY = meth$manifest[chr=='Y',index]
length(chrY)

detP = meth$detP[chrY,]
detP = colSums(detP<0.01,na.rm=TRUE)

pdf(file=paste("Ewastools_SexChrom_",date,".pdf",sep=""))
boxplot(split(detP,pheno$A1_SEX),ylab="# of detected Y chromosome probes")
split(detP,pheno$A1_SEX) %>% sapply(median)
dev.off()

sexpred = check_sex(meth) 
PredSex = predict_sex(sexpred$X, sexpred$Y)

pheno$predSex_X = sexpred$X
pheno$predSex_Y = sexpred$Y
pheno$predSex = PredSex

table(pheno$A1_SEX,pheno$predSex)

pheno$predSex_conv[pheno$predSex=="m"]=1
pheno$predSex_conv[pheno$predSex=="f"]=2
pheno$exclude[pheno$A1_SEX!=pheno$predSex_conv] = TRUE

sum(pheno$exclude)

#create beta matrix
beta = dont_normalize(meth)

# -------------------------
#SNP outliers
snps = meth$manifest[probe_type=="rs"]$index
genotypes = call_genotypes(beta[snps,],learn=FALSE)
pheno$SNPoutlier = snp_outliers(genotypes)

pdf(file=paste("Ewastools_snp_outlier_all_",date,".pdf",sep=""))
stripchart(pheno$SNPoutlier,method="jitter",pch=4)
abline(v=-4,lty="dotted",col=2)
dev.off()

pheno[pheno$SNPoutlier>=-4.0,.(SampleID,Sample.ID,Cohort)]

pheno[SNPoutlier > -4,exclude:=TRUE]
sum(pheno$exclude)

# ------------------------
#Check for duplicates
pheno$donor_id = enumerate_sample_donors(genotypes)

pheno[,n:=.N,by=donor_id]
pheno[n>1,.(SampleID,Sample.ID,Cohort,donor_id)]
pheno[n>1,exclude:=TRUE]

sum(pheno$exclude)

pheno = pheno[exclude==FALSE,]
dim(pheno)

write.csv(pheno, "EWastools_pheno.csv",row.names=FALSE)
