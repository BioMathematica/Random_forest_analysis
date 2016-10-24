#Make SG_otu_tax.RData

setwd("~/Dropbox (Second Genome)/CRC_wBaylor/Stool_Abundance_Table_SG/")
filedir <- "/Users//manasishah/Documents/Second_Genome/SG_rel_abun_RF"
otu_file <- read.csv("~/Dropbox (Second Genome)/CRC_wBaylor/Stool_Abundance_Table_SG/SG_CRC_FR_OTU_abundanceTable.csv", header=TRUE, row.names=1)
map_file <- "~/Dropbox (Second Genome)/CRC_wBaylor/Stool_Metadata_Table/Fecal_SG_MP.csv"
tax_file <- read.csv("~/Dropbox (Second Genome)/CRC_wBaylor/Stool_Abundance_Table_SG/SG_fecal9_otu_taxonomy.csv", header=TRUE)

# Check/remove quotes from this file.
library(phyloseq)
# importing the .csv otu_table-tax_table file.
SG_otu_tax_table <- otu_table(otu_file, taxa_are_rows = TRUE)
library("data.table")
#fread fast and friendly file finagler
map = fread(map_file, verbose=TRUE)
mapdf = as(map, "data.frame")

# Fix variable and sample names
colnames(mapdf)[1] <- "SampleID"
head(mapdf)
#Make syntactically valid names out of character vectors
sample_names(SG_otu_tax_table) <- make.names(sample_names(SG_otu_tax_table))
#Make the row names of mapdf == SampleID names/#Make syntactically valid names out of character vectors
rownames(mapdf) <- make.names(mapdf$SampleID)
mapdf <- sample_data(mapdf)
#Append the mapping file data as sample data for the phyloseq object
class(mapdf)
SG_otu_tax_table = merge_phyloseq(SG_otu_tax_table, mapdf)

#Load SG_phyloseq object with Taxonomy
load("~/SG_otu_tax.RData")
SG_otu_tax_CRC <- subset_samples(SG_otu_tax_CRC, !(Study == "Brim_V13_454"))
#Getting only CRC and control cases
SG_otu_tax_CRC <- subset_samples(SG_otu_tax, disease_stat%in%c("carcinoma", "control") )
SG_otu_tax_CRC <- subset_samples(SG_otu_tax_CRC, !(Treatment == "Yes"))
#Excluding controls from the study that had only adenoma and controls

#Total 339 samples retained, 5% ~ 17 samples
#Apply the prevalence filter, keep OTUs that occur in >= 17 samples #OTUs that occur in 5% of the samples
#n_to_keep = percent_to_keep*nsamples(SG_otu_tax_CRC)
#keepOTUs = apply(as(otu_table(SG_otu_tax_CRC), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= n_to_keep)
keepOTUs = apply(as(otu_table(SG_otu_tax_CRC), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= 17L)
SG_otu_tax_filt_5 <- prune_taxa(keepOTUs, SG_otu_tax_CRC)
#sgplot_filter_effect_library(SG_otu_tax_filt_5)
#sgplot_filter_effect_features(SG_otu_tax_filt_5)
SG_rel_filt <- transform_sample_counts(SG_otu_tax_filt_5, function(x) x/sum(x))
SG_filt_df <- data.frame(otu_table(SG_rel_filt))
SG_filt_t <- t(SG_filt_df)
SG_rel <- data.frame(SG_filt_t)
SG_filt5 <- cbind(SampleID = rownames(SG_rel), SG_rel)

#Get rowmeans for all OTU relative abundances to be used in dESeq2 figure, merge mean_rel_abun with top DESeq2 results

SG_filt_df$mean_rel_abun <- rowMeans(SG_filt_df)
SG_filt_df$OTU <- rownames(SG_filt_df)

df<- data.frame(sample_data(SG_otu_tax))
disease_stat <- subset(df, select=c("SampleID", "disease_stat", "Study"))
disease_stat <- subset(disease_stat, !(Study %in% "Brim_V13_454" ))
disease_CRC = disease_stat[disease_stat$disease_stat %in% c("carcinoma", "control"), ]
disease_CRC$disease_stat <- factor(disease_CRC$disease_stat)

SG_filt5 <- merge(SG_filt5, disease_CRC, by="SampleID")
row.names(SG_filt5) <- SG_filt5$SampleID
SG_filt5$SampleID <- NULL
SG_filt5$Study <- NULL


# ImpMeasure_filt5<-data.frame(varImp(SG_filt5_RF)$importance)
# Imp_SG_nzv <- data.frame(varImp(SG_rel_microbiome_only)$importance)
# ImpMeasure_filt5$OTU <- row.names(ImpMeasure_filt5)
# Imp_SG_nzv$OTU <- row.names(Imp_SG_nzv)

#commonOTUs <- merge(ImpMeasure_filt5, Imp_SG_nzv, by="OTU", all.x = F, all.y = F)

#commonOTUs <- order(commonOTUs$Overall.x)

############################################################
# 
# #Get SG OTU tax wit only CRC and control as disease status
# #Getting only CRC and control cases
# SG_otu_tax_CRC <- subset_samples(SG_otu_tax, disease_stat%in%c("carcinoma", "control") )
# #Excluding controls from the study that had only adenoma and controls
# SG_otu_tax_CRC <- subset_samples(SG_otu_tax_CRC, !(Study == "Brim_V13_454"))

# minusZack = subset_samples(SG_otu_tax, !(Study == "Zack_V4_MiSeq"))

studies = c("Zack_V4_MiSeq", "WuZhu_V3_454", "Wang_V3_454", "Chen_V13_454", "Zeller_V4_MiSeq", "Weir_V4_454")

studyList_minus = lapply(studies, function(x, input_physeq = SG_otu_tax_CRC){
  prune_samples(get_variable(input_physeq, "Study") != x, input_physeq)
})
names(studyList_minus) <- studies

#For single studies:
studyListsingle = lapply(studies, function(x, input_physeq = SG_otu_tax_CRC){
  prune_samples(get_variable(input_physeq, "Study") == x, input_physeq)
})

names(studyListsingle) <- studies

#Replace the  == sign in the above loop when you want to run a loop on single studies instead

prevalence_norm_fun = function(physeq, percent_to_keep){
  n_to_keep = percent_to_keep*nsamples(physeq)
  keep = apply(as(otu_table(physeq), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= n_to_keep)
  ps = prune_taxa(keep, physeq)
  ps = transform_sample_counts(ps, function(x) x/sum(x))
  otutab = t(otu_table(ps))
  otutab <- data.frame(otutab)
  otutab$disease_stat <- factor(get_variable(ps, "disease_stat"))
  return(otutab)
}

nMinus1_list = lapply(studyList_minus, prevalence_norm_fun, percent_to_keep = 0.05)
names(nMinus1_list) <- names(nMinus1_list)

nSingleStudy = lapply(studyListsingle, prevalence_norm_fun, percent_to_keep = 0.05)
names(nSingleStudy) <- names(nSingleStudy)

####################################################################################################################################################################################


#All RF's:

cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, savePred=T, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
SG_filt5_RF= train(disease_stat ~ ., SG_filt5, method = "rf", metric="ROC", trControl = cvCtrl, tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(SG_filt5)-1)))))

trainObject = SG_filt5_RF
tune.best = trainObject$bestTune
predMat = trainObject$pred
for(t in colnames(tune.best)){
  predMat = predMat[which(predMat[, t] == tune.best[, t]), ]
}
SG_filt5_conf <- confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")

cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, savePred=T, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
#newGrid = expand.grid(mtry = c((0.5, 1, 1.5, 1.75, 2)*(sqrt()))

resList_minus1 = lapply(nMinus1_list, function(x, trc){
  train(disease_stat ~. , x, method = "rf", metric = "ROC", trControl = trc, ntree = 1000, 
        tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(x)-1))))) 
        }, trc = cvCtrl)

confList_minus1 = lapply(resList_minus1, function(x){
  predMat = x$pred
  predMat = predMat[which(predMat[, "mtry"] == x$bestTune$mtry), ]
  confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")
})

####################################################################################################################################################################################

resList_single = lapply(nSingleStudy, function(x, trc){
  train(disease_stat ~. , x, method = "rf", metric = "ROC", trControl = trc, ntree = 1000, 
        tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2, 2.5, 3.0)*(sqrt(ncol(x)-1))))) 
}, trc = cvCtrl)

confList_single = lapply(resList_single, function(x){
  predMat = x$pred
  predMat = predMat[which(predMat[, "mtry"] == x$bestTune$mtry), ]
  confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")
})


####################################################################################################################################################################################
#Repeat Random Forest for n-1 leaving Weir out
SG_CRC_minusWeir <- subset_samples(SG_otu_tax_CRC, !(Study == "Weir_V4_454"))

studies_minusWeir = c("Zack_V4_MiSeq", "WuZhu_V3_454", "Wang_V3_454", "Chen_V13_454", "Zeller_V4_MiSeq")

studyList_minusWeir = lapply(studies_minusWeir, function(x, input_physeq = SG_CRC_minusWeir){
  prune_samples(get_variable(input_physeq, "Study") != x, input_physeq)
})
names(studyList_minusWeir) <- studies_minusWeir

prevalence_norm_fun = function(physeq, percent_to_keep){
  n_to_keep = percent_to_keep*nsamples(physeq)
  keep = apply(as(otu_table(physeq), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= n_to_keep)
  ps = prune_taxa(keep, physeq)
  ps = transform_sample_counts(ps, function(x) x/sum(x))
  otutab = t(otu_table(ps))
  otutab <- data.frame(otutab)
  otutab$disease_stat <- factor(get_variable(ps, "disease_stat"))
  return(otutab)
}

nMinusWeir1_list = lapply(studyList_minusWeir, prevalence_norm_fun, percent_to_keep = 0.05)
names(nMinusWeir1_list_list) <- names(nMinusWeir1_list)


cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, savePred=T, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
#newGrid = expand.grid(mtry = c((0.5, 1, 1.5, 1.75, 2)*(sqrt()))

resList_minusWeir1 = lapply(nMinusWeir1_list, function(x, trc){
  train(disease_stat ~. , x, method = "rf", metric = "ROC", trControl = trc, ntree = 1000, 
        tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(x)-1))))) 
}, trc = cvCtrl)

confList_minusWeir1 = lapply(resList_minusWeir1, function(x){
  predMat = x$pred
  predMat = predMat[which(predMat[, "mtry"] == x$bestTune$mtry), ]
  confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")
})


####################################################################################################################################################################################
#access the same value from the same function and loop it for the list to get desired values in a dataframe
####################################################################################################################################################################################

confList_single$Zack_V4_MiSeq$overall
confList_single$Zack_V4_MiSeq$positive
confList_single$Zack_V4_MiSeq$table
confList_single$Zack_V4_MiSeq$byClass
resList_single$Zack_V4_MiSeq$results$ROC

confList_minus1$Zack_V4_MiSeq$byClass
confList_minus1$Zeller_V4_MiSeq$byClass
confList_minus1$Weir_V4_454$byClass
confList_minus1$WuZhu_V3_454$byClass
confList_minus1$Chen_V13_454$byClass

#Get variable importance as a data frame

ImpMeasure_CRC_SG <-data.frame(varImp(SG_filt5_RF)$importance)
ImpMeasure_CRC_SG$OTU <- rownames(ImpMeasure_CRC_SG)
ImpMeasure_CRC_SG <- merge(ImpMeasure_CRC_SG, SG_tax_table[, c("OTU", "Phylum", "Family", "Genus", "Species", "Strain") ], by = "OTU")
ImpMeasure_CRC_SG <- ImpMeasure_CRC_SG[order(-ImpMeasure_CRC_SG$Overall),]

#Random forest for adenoma











