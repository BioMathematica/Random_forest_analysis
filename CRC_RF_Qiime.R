library(doMC)
ncores = detectCores()
registerDoMC(cores = ncores)
install.packages("caret")
install.packages("randomForest")
install.packages("e1071")
install.packages("pROC")
install.packages("glmnet")
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(glmnet)

#Load Q_phyloseq object with Taxonomy
load("~/Q_otu_tax.RData")
#Getting only CRC and control cases
Q_otu_tax_CRC <- subset_samples(Q_otu_tax, disease_stat%in%c("carcinoma", "control") )

#Excluding controls from the study that had only adenoma and controls
Q_otu_tax_CRC <- subset_samples(Q_otu_tax_CRC, !(Study == "Brim_V13_454"))
Q_otu_tax_CRC <- subset_samples(Q_otu_tax_CRC, !(Treatment == "Yes"))

#Total 424 samples retained, 5% ~ 21.8 samples
percent_to_keep <- 0.05
n_to_keep = percent_to_keep*nsamples(Q_otu_tax_CRC)
keepOTUs = apply(as(otu_table(Q_otu_tax_CRC), "matrix"), 1, function(x) sum(x >= 1L, na.rm = TRUE) >= n_to_keep)
Q_otu_tax_filt_5 <- prune_taxa(keepOTUs, Q_otu_tax_CRC)
#sgplot_filter_effect_library(Q_otu_tax_filt_5)
#sgplot_filter_effect_features(Q_otu_tax_filt_5)
Q_rel_filt <- transform_sample_counts(Q_otu_tax_filt_5, function(x) x/sum(x))
Q_filt_df <- data.frame(otu_table(Q_rel_filt))
Q_filt_t <- t(Q_filt_df)
Q_rel <- data.frame(Q_filt_t)
Q_filt5 <- cbind(SampleID = rownames(Q_rel), Q_rel)
Q_filt5_orig <- Q_filt5

#Get rowmeans for all OTU relative abundances to be used in dESeq2 figure, merge mean_rel_abun with top DESeq2 results

# Q_filt_df$mean_rel_abun <- rowMeans(Q_filt_df)
# Q_filt_df$OTU <- rownames(Q_filt_df)

df<- data.frame(sample_data(Q_otu_tax), stringsAsFactors = FALSE)
disease_stat <- subset(df, select=c("SampleID", "disease_stat", "Study"))
disease_stat <- subset(disease_stat, !(Study %in% "Brim_V13_454" ))
disease_CRC = disease_stat[disease_stat$disease_stat %in% c("carcinoma", "control"), ]
disease_CRC$disease_stat <- as.factor(disease_CRC$disease_stat)
Q_filt5 <- merge(Q_filt5, disease_CRC, by="SampleID", all.x = F, all.y = F)
Q_filt5$SampleID <- NULL
Q_filt5$Study <- NULL

####################################################################################################################################################################################

Q_micro_V_plat <- subset(df, select=c("SampleID", "disease_stat", "Study", "target_gene", "platform"))
Q_micro_V_plat <- subset(Q_micro_V_plat, !(Study %in% "Brim_V13_454" ))
Q_micro_V_plat = Q_micro_V_plat[Q_micro_V_plat$disease_stat %in% c("carcinoma", "control"), ]
Q_micro_V_plat[,c("target_gene","disease_stat")] = lapply(Q_micro_V_plat[, c("target_gene","disease_stat")], factor)
Q_micro_V_plat$disease_stat <- as.factor(Q_micro_V_plat$disease_stat)
Q_micro_V_plat$target_gene <- as.factor(Q_micro_V_plat$target_gene)
Q_micro_V_plat$platform <- factor(Q_micro_V_plat$platform)
Q_filt5_V_plat <- merge(Q_filt5_orig, Q_micro_V_plat, by="SampleID", all.x = F, all.y = F)
row.names(Q_filt5_V_plat) <- Q_filt5_V_plat$SampleID
Q_filt5_V_plat$SampleID <- NULL
Q_filt5_V_plat$Study <- NULL


################################################################################################################################################################################################################################################
#Trial model using near zero variance and highly correlated predictor

nzv <- nearZeroVar(Q_filt5_orig)
#973 features prior to filtering
Q_filt5_orig <- Q_filt5_orig[, -nzv]
dim(Q_filt5_orig)
#695 for 424 samples after nzv filtering
Q_filt5_orig$disease_stat <- NULL
row.names(Q_filt5_orig) <- Q_filt5_orig$SampleID
Q_filt5_orig$SampleID <- NULL
Q_filt5_orig<- as.matrix(Q_filt5_orig)
descrCor <-  cor(Q_filt5_orig)
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .95)
highCorr
#15 highly correlated variables 190 414 651 150 346 646 650 277 673   5 661
highlyCorasin <- findCorrelation(descrCor, cutoff = .95)
print(highlyCorasin)
filteredrel <- Q_filt5_orig[,-highlyCorasin]
Q_filt5_orig_abun <- data.frame(filteredrel)
Q_filt5_orig_abun <- cbind(SampleID = rownames(Q_filt5_orig_abun), Q_filt5_orig_abun)
Q_filt5_orig_abun2 <- Q_filt5_orig_abun
Q_filt5_orig_abun <- merge(Q_filt5_orig_abun, disease_CRC, by="SampleID")
row.names(Q_filt5_orig_abun) <- Q_filt5_orig_abun$SampleID
Q_filt5_orig_abun$SampleID <- NULL
Q_filt5_orig_abun$Study <- NULL
##########################################################################################################################################################


studies = c("Zack_V4_MiSeq", "WuZhu_V3_454", "Wang_V3_454", "Chen_V13_454", "Zeller_V4_MiSeq", "Weir_V4_454", "Pascual_V13_454", "Flemer_V34_MiSeq") 

studyList_minus_Q = lapply(studies, function(x, input_physeq = Q_otu_tax_CRC){
  prune_samples(get_variable(input_physeq, "Study") != x, input_physeq)
})
names(studyList_minus_Q) <- studies

#For single studies:
studyListsingle_Q = lapply(studies, function(x, input_physeq = Q_otu_tax_CRC){
  prune_samples(get_variable(input_physeq, "Study") == x, input_physeq)
})


names(studyListsingle_Q) <- studies

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

nMinus1_list_Q = lapply(studyList_minus_Q, prevalence_norm_fun, percent_to_keep = 0.05)
names(nMinus1_list_Q) <- names(nMinus1_list_Q)

nSingleStudy_Q = lapply(studyListsingle_Q, prevalence_norm_fun, percent_to_keep = 0.05)
names(nSingleStudy_Q) <- names(nSingleStudy_Q)

####################################################################################################################################################################################
nSingleStudy_Q[["Q_filt5"]] <- Q_filt5
nSingleStudy_Q[["Q_filt5_V_plat"]] <- Q_filt5_V_plat
nSingleStudy_Q[["Q_filt5_orig_abun"]] <- Q_filt5_orig_abun

#All RF's:

cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, savePred=T, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
Q_filt5_RF= train(disease_stat ~ ., Q_filt5, method = "rf", metric="ROC", trControl = cvCtrl, tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(Q_filt5)-1)))))

trainObject = Q_filt5_RF
tune.best = trainObject$bestTune
predMat = trainObject$pred
for(t in colnames(tune.best)){
  predMat = predMat[which(predMat[, t] == tune.best[, t]), ]
}
Q_filt5_conf <- confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")


#Q_filt5_VariableRegion_SequencingPlatform
Q_filt5_V_plat_RF= train(disease_stat ~ ., Q_filt5_V_plat, method = "rf", metric="ROC", trControl = cvCtrl, tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(Q_filt5)-1)))))

trainObject = Q_filt5_V_plat_RF
tune.best = trainObject$bestTune
predMat = trainObject$pred
for(t in colnames(tune.best)){
  predMat = predMat[which(predMat[, t] == tune.best[, t]), ]
}
Q_filt5_V_plat_conf <- confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")


Q_filt5_nzv_highCorr <- train(disease_stat ~ ., Q_filt5_orig_abun, method = "rf", metric="ROC", trControl = cvCtrl, tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(Q_filt5)-1)))))

trainObject = Q_filt5_nzv_highCorr
tune.best = trainObject$bestTune
predMat = trainObject$pred
for(t in colnames(tune.best)){
  predMat = predMat[which(predMat[, t] == tune.best[, t]), ]
}
Q_filt5_nzv_highCorr_conf <- confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")



cvCtrl = trainControl(method = "repeatedcv",number = 10, repeats = 5, savePred=T, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
resList_minus1_Q = lapply(nMinus1_list, function(x, trc){
  train(disease_stat ~. , x, method = "rf", metric = "ROC", trControl = trc, ntree = 1000, 
        tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2)*(sqrt(ncol(x)-1))))) 
}, trc = cvCtrl)

confList_minus1_Q = lapply(resList_minus1_Q, function(x){
  predMat = x$pred
  predMat = predMat[which(predMat[, "mtry"] == x$bestTune$mtry), ]
  confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")
})



####################################################################################################################################################################################

resList_single_Q = lapply(nSingleStudy_Q, function(x, trc){
  train(disease_stat ~. , x, method = "rf", metric = "ROC", trControl = trc, ntree = 1000, 
        tuneGrid = expand.grid(mtry = ceiling(c(0.5, 1, 1.5, 1.75, 2, 2.5, 3.0)*(sqrt(ncol(x)-1))))) 
}, trc = cvCtrl)

confList_single_Q = lapply(resList_single_Q, function(x){
  predMat = x$pred
  predMat = predMat[which(predMat[, "mtry"] == x$bestTune$mtry), ]
  confusionMatrix(predMat$pred, predMat$obs, positive = "carcinoma")
})



####################################################################################################################################################################################
#access the same value from the same function and loop it for the list to get desired values in a dataframe
####################################################################################################################################################################################

#Get variable importance as a data frame

ImpMeasure_CRC_Q_V_plat <-data.frame(varImp(Q_filt5_V_plat_RF)$importance)
ImpMeasure_CRC_Q_microbiome_only <-data.frame(varImp(Q_filt5_RF)$importance)
ImpMeasure_CRC_Q$OTU <- rownames(ImpMeasure_CRC_Q)
ImpMeasure_CRC_Q <- merge(ImpMeasure_CRC_Q, Q_tax[, c("OTU", "Phylum", "Family", "Genus", "Species", "Strain") ], by = "OTU")
ImpMeasure_CRC_Q <- ImpMeasure_CRC_Q[order(-ImpMeasure_CRC_Q$Overall),]

# Aggregate data results of the n-1 and single study random forest

conf_dt_Q = function(i, confList_minus1_Q){
  x = confList_minus1[[i]]$byClass
  dt1 = data.table(Type = names(x), Value = x)
  dt1[, Study := i]
  return(dt1)
}

conf_dt_single_Q = function(i, confList_single_Q){
  x = confList_single_Q[[i]]$byClass
  dt1 = data.table(Type = names(x), Value = x)
  dt1[, Study := i]
  return(dt1)
}

confdtminus1_Q = rbindlist(lapply(names(confList_minus1_Q), conf_dt_Q, confList_minus1_Q))
confdtsingle_Q = rbindlist(lapply(names(confList_single_Q), conf_dt_single_Q, confList_single_Q))


#Random forest for adenoma



