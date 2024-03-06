library(WGCNA)
library(flashClust)
library(readr)
library(datasets)
library(ggplot2)
library(tidyverse)

#####WGCNA Analysis#####
#Male mPFC
##Import normalized counts file that was created from DESeq2
datExpr_M <- read_csv("./DESeq2-mPFC-M-fullmodel.csv")
#An 'X' needs to be added to the top row of the first column of the file before reformatting the file
colnames(datExpr_M)[1] <- "X"
View(datExpr_M)

#Reformatting file for use in WGCNA
datExpr_M_matrix <- dat_Expr_M[ , -1]
row.names(datExpr_M_matrix) <- dat_Expr_M$X
datExpr_M_matrix <- as.matrix(datExpr_M_matrix)

datExpr_M_log <- log2(datExpr_M_matrix+1) #Log transform the count data

datExpr_M_log = as.data.frame(t(datExpr_M_log)) #Samples are rows and genes are columns
dim(datExpr_M_log)

#Run this to check if there are gene outliers (there always are)
gsg = goodSamplesGenes(datExpr_M_log, verbose = 3)
gsg $allOK

#Remove gene outliers (low counts and/or variability)

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr_M_log)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_M_log)[!gsg$goodSamples], collapse=", ")))
  datExpr_M_log= datExpr_M_log[gsg$goodSamples, gsg$goodGenes]
}

#Import file that contains the study variables
#Prenatal treatment group is dummy coded, all study variables need to be in numeric format otherwise the WGCNA code will not work
datTraits_M <- read_csv("./data/datTraits_mPFC_M.csv", 
                        col_types = cols(ID = col_character(), 
                                         Group_1 = col_number(), Group_2 = col_number(), 
                                         Group_3 = col_number(), Nest_attendance_P1_5_scaled = col_number(), 
                                         LG_P1_5_scaled = col_number()))
View(datTraits_M)

#Reformatting file
datTraits_M_matrix <- datTraits_M[ , -1]
row.names(datTraits_M_matrix) <- datTraits_M$ID

head(datTraits_M)

#form a data frame analogous to the count data
dim(datTraits_M)
table(rownames(datTraits_M_matrix)==rownames(datExpr_M_log)) #should return TRUE if datasets align correctly

#####Cluster count data#####

A = adjacency(t(datExpr_M_log),type="signed") #this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 #standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

#Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits_M_matrix,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits_M_matrix))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap") #no sample outliers were removed in this study

#Soft power was determined by a prior analysis with male mPFC data
enableWGCNAThreads()
softPower = 8
adjacency = adjacency(datExpr_M_log, power = softPower, type = "signed")

##Construct Networks

#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1-TOM


##Generate co-expressed eigengene modules

#Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#Minimum number of genes to cluster into a module was set to 30
minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)

#Create module eigengene lists
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr_M_log, colors= dynamicColors,softPower = 8)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")

#Plot dendrogram with module colors below
plotDendroAndColors(geneTree, cbind(dynamicColors), c("Eigengene Modules"), rowText=names(moduleLabels), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

#####Statistical Analyses#####

#Define number of genes and samples
nGenes = ncol(datExpr_M_log)
nSamples = nrow(datExpr_M_log)

##Import and reformat study traits file for statistical analysis
mPFC_design_M <- read_excel("./data/mPFC_tagseq_variables_M.xlsx")
View(mPFC_design_M)

#Factor prenatal BP groups, control group first
mPFC_design_M$Group <- factor(mPFC_design_M$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

#Scaling maternal care measures so 0 is approximately the mean and 1 is approximately the stdev
mPFC_design_M$Nest_attendance_P1_5_scaled <- (mPFC_design_M$Nest_attendance_P1_5-1600)/400
mPFC_design_M$LG_P1_5_scaled <- (mPFC_design_M$LG_P1_5-300)/100

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_M_log, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
colnames(MEs) #ME color names

#Incorporate study variables into MEs data
MEs$Group <- mPFC_design_M$Group
MEs$LG_P1_5_scaled <- mPFC_design_M$LG_P1_5_scaled
MEs$Nest_attendance_P1_5_scaled <- mPFC_design_M$Nest_attendance_P1_5_scaled
MEs$Group <- factor(MEs$Group, levels = c('Corn Oil', '50 μg/kg BPA', '50 μg/kg Mixed BP', '150 μg/kg Mixed BP'))

##Significance testing for prenatal BP treatment group, maternal care, and their interactions
#Outlier/insignificqnt modules were excluded
moduleTraitLm = lm(MElightcyan ~ Group + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEtan ~ Group + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEblack ~ Group + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEbrown ~ Group + LG_P1_5_scaled + Group:LG_P1_5_scaled, data = MEs)
moduleTraitLm = lm(MEyellow ~ Group + Nest_attendance_P1_5_scaled + Group:Nest_attendance_P1_5_scaled, data = MEs)

summary(moduleTraitLm)

#Exporting gene names within modules as a txt file
blackmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="black"])
yellowmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="yellow"])
lightcyanmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="lightcyan"])
tanmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="tan"])
brownmodule<-as.data.frame(names(datExpr_M_log)[moduleColors=="brown"])

write_delim(blackmodule, file = "M-mPFC-blackmodule.txt", delim = "\n")
write_delim(yellowmodule, file = "M-mPFC-yellowmodule.txt", delim = "\n")
write_delim(lightcyanmodule, file = "M-mPFC-lightcyan.txt", delim = "\n")
write_delim(tanmodule, file = "M-mPFC-tanmodule.txt", delim = "\n")
write_delim(brownmodule, file = "M-mPFC-brownmodule.txt", delim = "\n")

allgenes<-as.data.frame(names(datExpr_M_log)) #all genes analyzed for WGCNA
write_delim(allgenes, file = "M-mPFC-allgenes.txt", delim = "\n")

####Graphing MEs Data#####

#Subset the data by prenatal BP group so individual regression lines can be displayed
MEs_M_Control <- subset(MEs, Group == "Corn Oil")
MEs_M_BPA <- subset(MEs, Group == "50 μg/kg BPA")
MEs_M_BP_Low <- subset(MEs, Group == "50 μg/kg Mixed BP")
MEs_M_BP_High <- subset(MEs, Group == "150 μg/kg Mixed BP")

##ggplot2 scatterplot with points and regression lines colored by treatment group, Turquoise is used as an example module with Licking/grooming as the x variable
WGCNA_mPFC_M <- ggplot(MEs, aes(x = LG_P1_5_scaled, y = MEturquoise, color = Group)) +
  geom_point(size = 2.5,alpha = 0.55) + 
  scale_color_manual(values=c("#000000",'#739bd0',"#779F65",'#700064'), name = "Prenatal Treatment")+
  scale_y_continuous(name = "Turquoise Eigengene Value") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title = element_text(size = 16), legend.title = element_text(size=14), legend.text = element_text(size = 12))+
  scale_x_continuous(name = "Licking/Grooming Received per day from\nPND 1-5 (scaled)")+
  geom_smooth(data = MEs_M_Control, method = 'lm', se= FALSE, colour = "#000000", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BPA, method = 'lm', se= FALSE, colour = "#739bd0", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BP_Low, method = 'lm', se= FALSE, colour = "#779F65", fullrange = TRUE)+
  geom_smooth(data = MEs_M_BP_High, method = 'lm', se= FALSE, colour = "#700064", fullrange = TRUE)

WGCNA_mPFC_M
