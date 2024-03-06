library("DESeq2")
library(datasets)
library(ggplot2)
library(tidyverse)
library(readxl)

#####DEG Analysis####
#Female mPFC
#Import tagseq counts and study variable files to R
mPFC_tagseq_F <- read_excel("./data/mPFC_tagseq_counts_F.xlsx")
View(mPFC_tagseq_F)

mPFC_design_F <- read_excel("./data/mPFC_tagseq_variables_F.xlsx")
View(mPFC_design_F)

#Format file for DESeq2 analysis
mPFC_tagseq_F_matrix <- mPFC_tagseq_F[ , -1]
row.names(mPFC_tagseq_F_matrix) <- mPFC_tagseq_F$Gene_ID

#Factor prenatal BP groups, control group first
mPFC_design_F$Group <- factor(mPFC_design_F$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
mPFC_design_F$Nest_attendance_P1_5_scaled <- (mPFC_design_F$Nest_attendance_P1_5-1600)/400
mPFC_design_F$LG_P1_5_scaled <- (mPFC_design_F$LG_P1_5-300)/100

#Creating Deseq2 matrix using the full statistical model (prenatal BP treatment and postnatal maternal care)
ddsMatrix <- DESeqDataSetFromMatrix(
  countData = mPFC_tagseq_F_matrix, 
  colData = mPFC_design_F, 
  design = ~ 
    Group+
    Nest_attendance_P1_5_scaled+
    LG_P1_5_scaled+
    Group:Nest_attendance_P1_5_scaled+
    Group:LG_P1_5_scaled)

#Keep genes with 5 or more counts for 7 or more biological samples
keep <- rowSums(counts(ddsMatrix) >= 5) >= 7
ddsMatrix <- ddsMatrix[keep,]
dim(ddsMatrix)

#Estimate size factors, dispersion, normalize and perform negative binomial test
dds<-DESeq(ddsMatrix)

#Create file with normalized counts for use in WGCNA and/or plotting count data
ddsMatrix<- estimateSizeFactors(ddsMatrix)
normalized_counts <- counts(ddsMatrix, normalized=TRUE)
write.csv(normalized_counts, file="DESeq2-mPFC-F-fullmodel.csv", quote=FALSE)

###Collect results from the statistical testing, order by unadjusted pvalue
#Check the comparisons that were done
resultsNames(dds)

##Comparisons from full statistical model
#Main effects of prenatal BP group
res<-results(dds,name= "Group_50.μg.kg.BPA_vs_Corn.Oil")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group_50.μg.kg.Mixed.BP_vs_Corn.Oil")
resOrdered <- res[order(res_$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group_150.μg.kg.Mixed.BP_vs_Corn.Oil")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

#Main effects of maternal care
res<-results(dds,name= "Nest_attendance_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "LG_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

#Interactions between prenatal BP group and nest attendance
res<-results(dds,name= "Group50.μg.kg.BPA.Nest_attendance_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group50.μg.kg.Mixed.BP.Nest_attendance_P1_5_scaled")
resOrdered <- res[order(res_F$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group150.μg.kg.Mixed.BP.Nest_attendance_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

#Interactions between prenatal BP group and licking/grooming (LG)
res<-results(dds,name= "Group50.μg.kg.BPA.LG_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group50.μg.kg.Mixed.BP.LG_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

res<-results(dds,name= "Group150.μg.kg.Mixed.BP.LG_P1_5_scaled")
resOrdered <- res[order(res$pvalue),]
summary(res)
head(resOrdered) #Lists top 10 DEGs by unadjusted pvalue

#####Graphing interactions#####
plotExpr_F <- read_csv("./DESeq2-mPFC-F-fullmodel.csv")

#An 'X' needs to be added to the top row of the first column of the file before reformatting the file
colnames(plotExpr_F)[1] <- "X"
View(plotExpr_F)

#Reformatting file for ggplot2
plotExpr_F_matrix <- plotExpr_F[ , -1]
row.names(plotExpr_F_matrix) <- plotExpr_F$X
plotExpr_F = as.data.frame(t(plotExpr_F_matrix))

#Incorporating study variables into the normalized counts file
plotExpr_F$Group <- mPFC_design_F$Group
plotExpr_F$LG_P1_5_scaled <- mPFC_design_F$LG_P1_5_scaled
plotExpr_F$Nest_attendance_P1_5_scaled <- mPFC_design_F$Nest_attendance_P1_5_scaled

#Any gene that is graphed needs to be in the numeric form, Esr1 is used as an example
plotExpr_F$`Esr1` <- as.numeric(plotExpr_F$`Esr1`)

#Subset the data by prenatal BP group so individual regression lines can be displayed
plotExpr_F_Control <- subset(plotExpr_F, Group == "Corn Oil")
plotExpr_F_BPA <- subset(plotExpr_F, Group == "50 μg/kg BPA")
plotExpr_F_BP_Low <- subset(plotExpr_F, Group == "50 μg/kg Mixed BP")
plotExpr_F_BP_High <- subset(plotExpr_F, Group == "150 μg/kg Mixed BP")

#Factor prenatal BP treatment groups, control group first
plotExpr_F$Group <- factor(plotExpr_F$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#ggplot2 scatterplot with points and regression lines colored by treatment group, Esr1 is used as an example gene with Nest attendance as the x variable
MBxBP_mPFC_F <- ggplot(plotExpr_F, aes(x = Nest_attendance_P1_5_scaled, y = Esr1, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#739bd0',"#779F65",'#700064'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Normalized Esr1 Counts") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Time Spent on Nest per day from\nPND 1-5 (scaled)")+
  geom_smooth(data = plotExpr_F_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = plotExpr_F_BPA, method = 'lm', se= FALSE, colour = "#739bd0", fullrange = TRUE)+
  geom_smooth(data = plotExpr_F_BP_Low, method = 'lm', se= FALSE, colour = "#779F65", fullrange = TRUE)+
  geom_smooth(data = plotExpr_F_BP_High, method = 'lm', se= FALSE, colour = "#700064", fullrange = TRUE)

MBxBP_mPFC_F