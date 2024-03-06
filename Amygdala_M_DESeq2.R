library("DESeq2")
library(datasets)
library(ggplot2)
library(tidyverse)
library(readxl)

#####DEG Analysis####
#Male Amygdala
#Import tagseq counts and study variable files to R
Amy_tagseq_M <- read_excel("./data/Amygdala_tagseq_counts_M.xlsx")
View(Amy_tagseq_M)

Amy_design_M <- read_excel("./data/Amygdala_tagseq_variables_M.xlsx")
View(Amy_design_M)

#Format file for DESeq2 analysis
Amy_tagseq_M_matrix <- Amy_tagseq_M[ , -1]
row.names(Amy_tagseq_M_matrix) <- Amy_tagseq_M$Gene_ID

#Factor prenatal BP groups, control group first
Amy_design_M$Group <- factor(Amy_design_M$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
Amy_design_M$Nest_attendance_P1_5_scaled <- (Amy_design_M$Nest_attendance_P1_5-1600)/400
Amy_design_M$LG_P1_5_scaled <- (Amy_design_M$LG_P1_5-300)/100

#Creating Deseq2 matrix using the full statistical model (prenatal BP treatment and postnatal maternal care)
ddsMatrix <- DESeqDataSetFromMatrix(
  countData = Amy_tagseq_M_matrix, 
  colData = Amy_design_M, 
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
write.csv(normalized_counts, file="DESeq2-Amygdala-M-fullmodel.csv", quote=FALSE)

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
plotExpr_M <- read_csv("./DESeq2-mPFC-M-fullmodel.csv")

#An 'X' needs to be added to the top row of the first column of the file before reformatting the file
colnames(plotExpr_M)[1] <- "X"
View(plotExpr_M)

#Reformatting file for ggplot2
plotExpr_M_matrix <- plotExpr_M[ , -1]
row.names(plotExpr_M_matrix) <- plotExpr_M$X
plotExpr_M = as.data.frame(t(plotExpr_M_matrix))

#Incorporating study variables into the normalized counts file
plotExpr_M$Group <- Amy_design_M$Group
plotExpr_M$LG_P1_5_scaled <- Amy_design_M$LG_P1_5_scaled
plotExpr_M$Nest_attendance_P1_5_scaled <- Amy_design_M$Nest_attendance_P1_5_scaled

#Any gene that is graphed needs to be in the numeric form, Esr1 is used as an example
plotExpr_M$`Esr1` <- as.numeric(plotExpr_M$`Esr1`)

#Subset the data by prenatal BP group so individual regression lines can be displayed
plotExpr_M_Control <- subset(plotExpr_M, Group == "Corn Oil")
plotExpr_M_BPA <- subset(plotExpr_M, Group == "50 μg/kg BPA")
plotExpr_M_BP_Low <- subset(plotExpr_M, Group == "50 μg/kg Mixed BP")
plotExpr_M_BP_High <- subset(plotExpr_M, Group == "150 μg/kg Mixed BP")

#Factor prenatal BP treatment groups, control group first
plotExpr_M$Group <- factor(plotExpr_M$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#ggplot2 scatterplot with points and regression lines colored by treatment group, Esr1 is used as an example gene with Nest attendance as the x variable
MBxBP_Amy_M <- ggplot(plotExpr_M, aes(x = Nest_attendance_P1_5_scaled, y = Esr1, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#739bd0',"#779F65",'#700064'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Normalized Esr1 Counts") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Time Spent on Nest per day from\nPND 1-5 (scaled)")+
  geom_smooth(data = plotExpr_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = plotExpr_M_BPA, method = 'lm', se= FALSE, colour = "#739bd0", fullrange = TRUE)+
  geom_smooth(data = plotExpr_M_BP_Low, method = 'lm', se= FALSE, colour = "#779F65", fullrange = TRUE)+
  geom_smooth(data = plotExpr_M_BP_High, method = 'lm', se= FALSE, colour = "#700064", fullrange = TRUE)

MBxBP_Amy_M