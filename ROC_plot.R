library(caret)
#load("css_results.RData")
resamp <- resamples(list(Clinical_RF_css, SG_rel_abun_microbiome, SG_rel_both), modelNames = c("Clinical Only", "Strain Only", "Clinical + Microbiome"))
resamp.sum <- summary(resamp)
write.csv(resamp.sum$statistics, "RF_resultMetrics.csv")
library(lattice)
bwplot(resamp)

#plot roc curves

library(pROC)
predMat <- both_rf_css$pred
predMat <- predMat[predMat$mtry == both_rf_css$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
pred = predMat$carcinoma
bothROC = roc(response=resp, predictor=pred)

predMat <- strain_RF_css$pred
predMat <- predMat[predMat$mtry == both_rf_css$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
pred = predMat$carcinoma
strainROC = roc(response=resp, predictor=pred)

predMat <- Clinical_RF_css$pred
predMat <- predMat[predMat$mtry == Clinical_RF_css$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
pred = predMat$carcinoma
clinicalROC = roc(response=resp, predictor=pred)

predMat <- Q_clin_micro_train$pred
predMat <- predMat[predMat$mtry == Q_clin_micro_train$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
pred = predMat$carcinoma
Q_clin_micro_roc = roc(response=resp, predictor=pred)


plot.roc(strainROC, col = "green")
plot.roc(bothROC, add=T, col = "tomato")
plot.roc(clinicalROC, add=T, col = "cornflowerblue")
legend("bottomright", legend = c("Clinical + Microbiome", "Strain Only", "Clinical Only"), col = c("tomato","green", "cornflowerblue"), bty = "n", pch = "-", cex=1, pt.cex = 1.5)


#********************************************************************************************************************************************

#nMinus1_list[7] <- list(SG_filt5_RF)

nMinus1_resamp <- resamples(resList_minus1, modelNames = c("Minus Zack_V4_MiSeq", "Minus WuZhu_V3_454", "Minus Wang_V3_454", "Minus Chen_V13_454", "Minus Zeller_V4_MiSeq", "Minus Weir_V4_454", "Minus Pascual_V13_454", "Minus Flemer_V34_MiSeq"))

SingleStudy_resamp <- resamples(resList_single["Zack_V4_MiSeq", "WuZhu_V3_454", "Wang_V3_454", "Chen_V13_454", "Zeller_V4_MiSeq", "Weir_V4_454", "Pascual_V13_454", "Flemer_V34_MiSeq"], modelNames = c("Only Zack_V4_MiSeq", "Only WuZhu_V3_454", "Only Wang_V3_454", "Only Chen_V13_454", "Only Zeller_V4_MiSeq", "Only Weir_V4_454", "Only Pascual_V13_454", "Only Flemer_V34_MiSeq"))

nMinus1_resamp_sum <- summary(nMinus1_resamp)
write.csv(nMinus1_resamp_sum$statistics, "n-1_Results.csv")
library(lattice)
bwplot(nMinus1_resamp)

#********************************************************************************************************************************************

ROC_list = lapply(resList_minus1, function(x){
  predMat = x$pred
  predMat <- predMat[predMat$mtry == x$bestTune[[1]], ]
  resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
  pred = predMat$carcinoma
  ROC = roc(response=resp, predictor=pred, plot = T)
  return(ROC)
})


ROC_list_single = lapply(resList_single, function(x){
  predMat = x$pred
  predMat <- predMat[predMat$mtry == x$bestTune[[1]], ]
  resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
  pred = predMat$carcinoma
  ROC = roc(response=resp, predictor=pred, plot = T)
  return(ROC)
})


#Plot the ROC curve
tiff("Fig3A.tiff", width = 11.5, height = 6, units = "in", res = 300,
     compression = "lzw", colortype = "true", family="Arial")
#theme_set(theme_bw(base_size = 10))
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
title(outer=outer,adj=0,main = list("3A", cex=1.1,col="black", 
                                   font=2)) 
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(ROC_list_basic$`RF_list_res$Clinical_ML`, col = "black", add=T)
plot.roc(ROC_list_basic$`resList_single$SG_filt5`, col = "chartreuse4", add=T)
plot.roc(ROC_list_basic$`resList_single_Q$Q_filt5`, col = "red", add=T)
plot.roc(ROC_list_basic$`RF_list_res$SG_clin_micro`, col = "blue", add=T)
plot.roc(ROC_list_basic$Q_clin_micro_train, col = "orange", add=T)
legend("bottomright", legend = expression("Clinical features only (n=156, AUC=79.7)", "Microbial features SS-UP (n=424, AUC=80.4)*", 
                                          "Microbial features QIIME-CR (n=424, AUC=76.6, DeLong's p=0.008)*", "Microbial+Clinical features SS-UP (n=156, AUC=91.3)^", 
                                          "Microbial+Clinical features QIIME-CR (n=156, AUC=82.4, DeLong's p=1.2*" ~ 10^-6 ~ ")^"), 
       col = c("black", "chartreuse4", "red", "blue", "orange"), bty = "n", pch = "---", cex=0.8, pt.cex = 1.75)
invisible(dev.off())

tiff("Fig3B.tiff", width = 11, height = 6, units = "in", res = 300,
     compression = "lzw", colortype = "true", family="Arial")
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
title(outer=outer,adj=0,main = list("3B", cex=1.1,col="black", 
                                    font=2)) 
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(ROC_list$Zack_V4_MiSeq, col = "#9F832D", add=T)
plot.roc(ROC_list$WuZhu_V3_454, add=T, col="#F1BA2F")
plot.roc(ROC_list$Wang_V3_454, add=T, col="#DF8236")
plot.roc(ROC_list$Chen_V13_454, add=T, col= "#01B5BB")
plot.roc(ROC_list$Zeller_V4_MiSeq, add=T, col="#154db3")
plot.roc(ROC_list$Weir_V4_454, add=T, col="#228b22")
plot.roc(ROC_list$Pascual_V13_454, add=T, col="#8F389E")
plot.roc(ROC_list$Flemer_V34_MiSeq, add=T, col="#E50E63")
plot.roc(ROC_list_single$SG_filt5, add=T, col="black")
legend("bottomright", legend = c("Minus Zackular_V4_MiSeq", "Minus WuZhu_V3_454 AUC=83.9, DeLong's p=0.005*", "Minus Wang_V3_454 AUC=75.8, DeLong's p=0.003^", "Minus Chen_V13_454", "Minus Zeller_V4_MiSeq", "Minus Weir_V4_454", "Minus Pascual_V13_454", "Minus Flemer_V34_MiSeq", "Total cohort SS-UP AUC=80.3*^"), col = c("#9F832D","#F1BA2F", "#DF8236", "#01B5BB", "#154db3", "#228b22", "#8F389E", "#E50E63","black"), bty = "n", pch = "---", cex=0.8, pt.cex = 2)

invisible(dev.off())

tiff("Fig3C.tiff", width = 9.9, height = 6, units = "in", res = 300,
     compression = "lzw", colortype = "true", family="Arial")
par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
title(outer=outer,adj=0,main = list("3C", cex=1.1,col="black", 
                                    font=2))
#plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(ROC_list_single$Zack_V4_MiSeq, col = "#9F832D", add=T)
plot.roc(ROC_list_single$WuZhu_V3_454, add=T, col="#F1BA2F")
plot.roc(ROC_list_single$Wang_V3_454, add=T, col="#DF8236")
plot.roc(ROC_list_single$Chen_V13_454, add=T, col="#01B5BB")
plot.roc(ROC_list_single$Zeller_V4_MiSeq, add=T, col="#154db3")
plot.roc(ROC_list_single$Weir_V4_454, add=T, col="#228b22")
plot.roc(ROC_list_single$Pascual_V13_454, add=T, col="#8F389E")
plot.roc(ROC_list_single$Flemer_V34_MiSeq, add=T, col="#E50E63")
plot.roc(ROC_list_single$SG_filt5, add=T, col="black")
legend("bottomright", legend = c("Zackular_V4_MiSeq", "WuZhu_V3_454", " Wang_V3_454", " Chen_V13_454", " Zeller_V4_MiSeq", " Weir_V4_454", " Pascual_V13_454", " Flemer_V34_MiSeq", "Total cohort SS-UP"), col = c("#9F832D", "#F1BA2F", "#DF8236", "#01B5BB", "#154db3", "#228b22", "#8F389E", "#E50E63", "black"), bty = "n", pch = "---", cex=0.8, pt.cex = 1.75)
invisible(dev.off())

save(ROC_list_single, ROC_list_basic, ROC_list, resList_minus1, resList_single, file="Figure3.RData")




roc.test(ROC_list_single$SG_filt5, ROC_list_single$WuZhu_V3_454)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Zeller_V4_MiSeq)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Zack_V4_MiSeq)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Chen_V13_454)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Weir_V4_454)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Pascual_V13_454)
roc.test(ROC_list_single$SG_filt5, ROC_list_single$Flemer_V34_MiSeq)


#ROC list - Clinical ML, SS-UP microbial ML, QIIME-CR microbial ML, Clin+micro SS-UP and QIIME CR

name_list_objects <- function(...){
  names <- as.list(substitute(list(...)))[-1L]
  result <- list(...)
  names(result) <- names
  result
}

plot_list_ROC_basic <- name_list_objects(resList_single$SG_filt5, resList_single_Q$Q_filt5, Q_clin_micro_train, RF_list_res$Clinical_ML, RF_list_res$SG_clin_micro)
#Name lists objects the same as the data frames

ROC_list_basic = lapply(plot_list_ROC_basic, function(x){
  predMat = x$pred
  predMat <- predMat[predMat$mtry == x$bestTune[[1]], ]
  resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
  pred = predMat$carcinoma
  ROC = roc(response=resp, predictor=pred, plot = T)
  return(ROC)
})

roc.test(ROC_list_basic$`resList_single$SG_filt5`, ROC_list_basic$`resList_single_Q$Q_filt5`) 
#p=0.008
roc.test(ROC_list_basic$`RF_list_res$Clinical_ML`, ROC_list_basic$`resList_single$SG_filt5`)
#p=0.76
roc.test(ROC_list_basic$Q_clin_micro_train, ROC_list_basic$`RF_list_res$SG_clin_micro`)
#p= 1.2*10-6





########################################################################################################################################

#QIIME n-1 and single study plots

ROC_list_Q = lapply(resList_minus1_Q, function(x){
  predMat = x$pred
  predMat <- predMat[predMat$mtry == x$bestTune[[1]], ]
  resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
  pred = predMat$carcinoma
  ROC = roc(response=resp, predictor=pred, plot = T)
  return(ROC)
})


ROC_list_single_Q = lapply(resList_single_Q, function(x){
  predMat = x$pred
  predMat <- predMat[predMat$mtry == x$bestTune[[1]], ]
  resp = factor(predMat$obs, ordered = T, levels = c("control", "carcinoma"))
  pred = predMat$carcinoma
  ROC = roc(response=resp, predictor=pred, plot = T)
  return(ROC)
})

par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
title(outer=outer,adj=0,main = list("A", cex=1.1,col="black", 
                                    font=2)) 

mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(ROC_list_Q$Zack_V4_MiSeq, col = "cyan4", add=T)
plot.roc(ROC_list_Q$WuZhu_V3_454, add=T, col="tomato")
plot.roc(ROC_list_Q$Wang_V3_454, add=T, col="purple4")
plot.roc(ROC_list_Q$Chen_V13_454, add=T, col="violet")
plot.roc(ROC_list_Q$Zeller_V4_MiSeq, add=T, col="#840000")
plot.roc(ROC_list_Q$Weir_V4_454, add=T, col="orange")
plot.roc(ROC_list_Q$Pascual_V13_454, add=T, col="chartreuse4")
plot.roc(ROC_list_Q$Flemer_V34_MiSeq, add=T, col="cornflowerblue")
plot.roc(ROC_list_single_Q$Q_filt5, add=T, col="black")
legend("bottomright", legend = c("Minus Zackular_V4_MiSeq", "Minus WuZhu_V3_454", "Minus Wang_V3_454 AUC=0.71, DeLong's p=0.0003", "Minus Chen_V13_454", "Minus Zeller_V4_MiSeq", "Minus Weir_V4_454", "Minus Pascual_V13_454", "Minus Flemer_V34_MiSeq", "Total cohort QIIME-CR AUC=0.77"), col = c("cyan4","tomato", "purple4", "violet", "#840000", "orange", "chartreuse4", "cornflowerblue", "black"), bty = "n", pch = "---", cex=0.8, pt.cex = 2)
#title(main = "QIIME-CR leave one study out random forest", cex.main=1, font.main=2)

roc.test(ROC_list_single$SG_filt5, ROC_list_single$WuZhu_V3_454)


par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
title(outer=outer,adj=0,main = list("B", cex=1.1,col="black", 
                                    font=2)) 
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(ROC_list_single_Q$Zack_V4_MiSeq, col = "blue", add=T)
plot.roc(ROC_list_single_Q$WuZhu_V3_454, add=T, col="tomato")
plot.roc(ROC_list_single_Q$Wang_V3_454, add=T, col="purple4")
plot.roc(ROC_list_single_Q$Chen_V13_454, add=T, col="violet")
plot.roc(ROC_list_single_Q$Zeller_V4_MiSeq, add=T, col="#840000")
plot.roc(ROC_list_single_Q$Weir_V4_454, add=T, col="orange")
plot.roc(ROC_list_single_Q$Pascual_V13_454, add=T, col="chartreuse4")
plot.roc(ROC_list_single_Q$Flemer_V34_MiSeq, add=T, col="cornflowerblue")
plot.roc(ROC_list_single_Q$Q_filt5, add=T, col="black")
legend("bottomright", legend = c("Zackular_V4_MiSeq", " WuZhu_V3_454", " Wang_V3_454", " Chen_V13_454", " Zeller_V4_MiSeq", " Weir_V4_454", " Pascual_V13_454", " Flemer_V34_MiSeq", "Total cohort QIIME-CR"), col = c("blue","tomato", "purple4", "violet", "#840000", "orange", "chartreuse4","cornflowerblue", "black"), bty = "n", pch = "---", cex=0.8, pt.cex = 1.75)
title(main = "QIIME-CR individual study", cex.main=1, font.main=2)





roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$WuZhu_V3_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Zeller_V4_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Zack_V4_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Chen_V13_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Weir_V4_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Pascual_V13_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Flemer_V34_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_Q$Wang_V3_454)

roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$WuZhu_V3_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Zeller_V4_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Zack_V4_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Chen_V13_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Weir_V4_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Pascual_V13_454)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Flemer_V34_MiSeq)
roc.test(ROC_list_single_Q$Q_filt5, ROC_list_single_Q$Wang_V3_454)



#######################################################################################################################################

#Variable Importance curves

Imp_SG_micro <- data.frame(varImp(resList_single$SG_filt5)$importance)
Imp_SG_micro$OTU <- rownames(Imp_SG_micro)
Imp_SG_micro <- merge(Imp_SG_micro, SG_tax_table[c("OTU", "Taxonomy")], by="OTU")
Imp_SG_micro <- Imp_SG_micro[order(Imp_SG_micro$Overall),] 
Imp_SG_micro_filt <- subset(Imp_SG_micro, Imp_SG_micro$Overall > 10)

Imp_SG_micro_filt2 <- Imp_SG_micro_filt

Imp_SG_micro_filt2$Taxonomy <-factor(Imp_SG_micro_filt$Taxonomy, levels=Imp_SG_micro_filt[order(Imp_SG_micro_filt$Overall), "Taxonomy"])

p1 <- ggplot(Imp_SG_micro_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw()  + scale_x_reverse() +ggtitle("Microbial feature ranking (n=344), SS-UP") + xlab("Overall feature importance")


Imp_SG_micro_clin <- data.frame(varImp(RF_list_res$SG_clin_micro)$importance)
Imp_SG_micro_clin$OTU <- rownames(Imp_SG_micro_clin)
Imp_SG_micro_clin <- merge(Imp_SG_micro_clin, SG_tax_table[c("OTU", "Taxonomy")], by="OTU", all.x = T)
Imp_SG_micro_clin <- Imp_SG_micro_clin[order(Imp_SG_micro_clin$Overall),] 
Imp_SG_micro_clin_filt <- subset(Imp_SG_micro_clin, Imp_SG_micro_clin$Overall > 10)

Imp_SG_micro_clin_filt <- Imp_SG_micro_clin_filt[order(-Imp_SG_micro_clin_filt$Overall), ]

Imp_SG_micro_clin_filt[2, 3] <- "FOBT_P"
Imp_SG_micro_clin_filt[4, 3] <- "FOBT_N"


Imp_SG_micro_clin_filt2 <- Imp_SG_micro_clin_filt

Imp_SG_micro_clin_filt2$Taxonomy <-factor(Imp_SG_micro_clin_filt$Taxonomy, levels=Imp_SG_micro_clin_filt[order(Imp_SG_micro_clin_filt$Overall), "Taxonomy"])

# Imp_SG_micro_clin_filt2 <- Imp_SG_micro_clin_filt2[order(-Imp_SG_micro_clin_filt2$Overall), ]
# 
# Imp_SG_micro_clin_filt2[2, 3] <- "FOBT_P"
# Imp_SG_micro_clin_filt2[4, 3] <- "FOBT_N"

p2 <- ggplot(Imp_SG_micro_clin_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw() + scale_x_reverse() + ggtitle("Clinical + Microbial feature ranking (n=156), SS-UP") + xlab("Overall feature importance")



Imp_SG_micro_minusWu <- data.frame(varImp(resList_minus1$WuZhu_V3_454)$importance)
Imp_SG_micro_minusWu$OTU <- rownames(Imp_SG_micro_minusWu)
Imp_SG_micro_minusWu <- merge(Imp_SG_micro_minusWu, SG_tax_table[c("OTU", "Taxonomy")], by="OTU", all.x = T)
Imp_SG_micro_minusWu <- Imp_SG_micro_minusWu[order(Imp_SG_micro_minusWu$Overall),] 
Imp_SG_micro_minusWu_filt <- subset(Imp_SG_micro_minusWu, Imp_SG_micro_minusWu$Overall > 10)

Imp_SG_micro_minusWu_filt2 <- Imp_SG_micro_minusWu_filt

Imp_SG_micro_minusWu_filt2$Taxonomy <-factor(Imp_SG_micro_minusWu_filt$Taxonomy, levels=Imp_SG_micro_minusWu_filt[order(Imp_SG_micro_minusWu_filt$Overall), "Taxonomy"])

p3 <- ggplot(Imp_SG_micro_minusWu_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw()  + scale_x_reverse() +ggtitle("Microbial cohort minus Wu_V13_454 feature ranking, SS-UP") + xlab("Overall feature importance")

#######################################################################################################################################

Q_tax_table <- as.data.frame(tax_table(Q_otu_tax))
Q_tax_table <- within(Q_tax_table, Taxonomy <- paste(Phylum,Family,Genus,Species, sep=';'))
Q_tax_table$OTU <- rownames(Q_tax_table)

Imp_Q_micro <- data.frame(varImp(resList_single_Q$Q_filt5)$importance)
Imp_Q_micro$OTU <- rownames(Imp_Q_micro)
Imp_Q_micro <- merge(Imp_Q_micro, Q_tax_table[c("OTU", "Taxonomy")], by="OTU")
Imp_Q_micro <- Imp_Q_micro[order(Imp_Q_micro$Overall),] 
Imp_Q_micro_filt <- subset(Imp_Q_micro, Imp_Q_micro$Overall > 10)

Imp_Q_micro_filt2 <- Imp_Q_micro_filt

Imp_Q_micro_filt2$Taxonomy <-factor(Imp_Q_micro_filt$Taxonomy, levels=Imp_Q_micro_filt[order(Imp_Q_micro_filt$Overall), "Taxonomy"])

p4 <- ggplot(Imp_Q_micro_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw()  + scale_x_reverse() +ggtitle("Microbial feature ranking (n=344), QIIME-CR") + xlab("Overall feature importance")

#######################################################################################################################################
library("gridExtra")
#ggsave(fsDESeq_plot_Control_Vs_adenoma, plot = "manasisfavoriteplot.pdf", width = 8, height = 10)
varImp_plot_list = list(p1, p2, p3, p4)
do.call(grid.arrange,  varImp_plot_list)
# ggplot(Imp_SG_micro, aes(x=Overall, y=Taxonomy)) + geom_point() + scale_x_reverse() + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black"))
# ggplot(Imp_SG_micro_filt, aes(x=reorder(Overall), y=Taxonomy)) + geom_point() + scale_x_reverse() + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black"))
# ggplot(Imp_SG_micro_filt, aes(x=(reorder(Overall)), y=Taxonomy)) + geom_point() + scale_x_reverse() + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black"))
# ggplot(Imp_SG, aes(x=(reorder(Overall)), y=Taxonomy)) + geom_point() + scale_x_reverse() + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black"))

#######################################################################################################################################

library(pROC)
predMat <- Q_paired_fb$pred
predMat <- predMat[predMat$mtry == Q_paired_fb$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("fecal_carcinoma", "biopsy_carcinoma"))
pred = predMat$biopsy_carcinoma
FB_ROC = roc(response=resp, predictor=pred)

library(pROC)
predMat <- Q_paired_biopsy$pred
predMat <- predMat[predMat$mtry == Q_paired_biopsy$bestTune[[1]], ]
resp = factor(predMat$obs, ordered = T, levels = c("carcinoma_adjacent", "biopsy_carcinoma"))
pred = predMat$biopsy_carcinoma
BB_ROC = roc(response=resp, predictor=pred)

par(mar=c(4,4,1,1))
plot(c(1,0),c(0,1), type='l', lty=3, xlim=c(1.01,0), ylim=c(-0.01,1.01), xaxs='i', yaxs='i', ylab='', xlab='')
mtext(side=2, text="Sensitivity", line=2.5)
mtext(side=1, text="Specificity", line=2.5)
plot.roc(FB_ROC, col = "cyan4", add=T)
plot.roc(BB_ROC, add=T, col="tomato")
legend("bottomright", legend = c("Fecal Vs Biopsy CRC Samples, AUC = 71.1", "CRC biopsy Vs paired tumor adjacent biopsies, AUC = 66.7"), col = c("cyan4","tomato"), bty = "n", pch = "---", cex=0.8, pt.cex = 1.75)
title(main = "Biopsy classifier", cex.main=1, font.main=2)


#Variable importance for the two models 
Q_tax_table <- as.data.frame(tax_table(Q_otu_tax))
Q_tax_table$OTU <- rownames(Q_tax_table)
Q_tax_table$Genus <- gsub("g__", "", Q_tax_table$Genus)
Q_tax_table$Genus <- gsub("^$", "unc", Q_tax_table$Genus)
Q_tax_table$Phylum <- gsub("p__", "", Q_tax_table$Phylum)
Q_tax_table$Phylum <- gsub("^$", "unc", Q_tax_table$Phylum)
Q_tax_table$Species <- gsub("s__", "", Q_tax_table$Species)
Q_tax_table$Species <- gsub("^$", "unc", Q_tax_table$Species)
Q_tax_table$Family <- gsub("f__", "", Q_tax_table$Species)
Q_tax_table$Family <- gsub("^$", "unc", Q_tax_table$Species)

Q_tax_table$Species[Q_tax_table$Species == "unc"] <- as.character(Q_tax_table$OTU[Q_tax_table$Species == "unc"])

#Q_tax_table$Species[Q_tax_table$Species == "unc"] <- Q_tax_table$OTU

Q_tax_table <- within(Q_tax_table, Taxonomy <- paste(Phylum,Genus,Species, sep=';'))

bb_imp <- data.frame(varImp(Q_paired_biopsy)$importance)
bb_imp$OTU <- rownames(bb_imp)
bb_imp$OTU <- gsub("X","OTU", bb_imp$OTU)
bb_imp <- merge(bb_imp, Q_tax_table[c("OTU", "Taxonomy")], by="OTU")
bb_imp <- bb_imp[order(-bb_imp$Overall),] 
bb_imp_filt <- subset(bb_imp, bb_imp$Overall > 10)

bb_imp_filt2 <- bb_imp_filt

bb_imp_filt2$Taxonomy <-factor(bb_imp_filt$Taxonomy, levels=bb_imp_filt[order(bb_imp_filt$Overall), "Taxonomy"])



bb_imp_filt2 <- bb_imp_filt

bb_imp_filt2$Taxonomy <-factor(bb_imp_filt$Taxonomy, levels=bb_imp_filt[order(bb_imp_filt$Overall), "Taxonomy"])

bb <- ggplot(head(bb_imp_filt2, 20), aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_set(theme_bw(base_size = 15))   + scale_x_reverse()  + xlab("Overall feature importance")


fb_imp <- data.frame(varImp(Q_paired_fb)$importance)
fb_imp$OTU <- rownames(fb_imp)
fb_imp$OTU <- gsub("X","OTU", fb_imp$OTU)
fb_imp <- merge(fb_imp, Q_tax_table[c("OTU", "Taxonomy")], by="OTU")
fb_imp <- fb_imp[order(-fb_imp$Overall),] 
fb_imp_filt <- subset(fb_imp, fb_imp$Overall > 10)

fb_imp_filt2 <- fb_imp_filt

fb_imp_filt2$Taxonomy <-factor(fb_imp_filt$Taxonomy, levels=fb_imp_filt[order(fb_imp_filt$Overall), "Taxonomy"])

fb <- ggplot(head(fb_imp_filt2, 20), aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_set(theme_bw(base_size = 15))   + scale_x_reverse()  + xlab("Overall feature importance")







FB_random_forest <- ggplot(Imp_Q_micro_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + 
geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw()  + 
scale_x_reverse() +ggtitle("Microbial feature ranking, fecal biopsy classifier") + xlab("Overall feature importance")

Imp_Q_micro <- data.frame(varImp(Q_paired_biopsy)$importance)
Imp_Q_micro$OTU <- rownames(Imp_Q_micro)
Imp_Q_micro$OTU <- gsub("X","OTU", Imp_Q_micro$OTU)
Imp_Q_micro <- merge(Imp_Q_micro, Q_tax_table[c("OTU", "Taxonomy")], by="OTU")
Imp_Q_micro <- Imp_Q_micro[order(Imp_Q_micro$Overall),] 
Imp_Q_micro_filt <- subset(Imp_Q_micro, Imp_Q_micro$Overall > 10)

Imp_Q_micro_filt2 <- Imp_Q_micro_filt

Imp_Q_micro_filt2$Taxonomy <-factor(Imp_Q_micro_filt$Taxonomy, levels=Imp_Q_micro_filt[order(Imp_Q_micro_filt$Overall), "Taxonomy"])
Imp_Q_micro_filt2 <- Imp_Q_micro_filt2[order(-Imp_Q_micro_filt2$Overall),] 
BB_random_forest <- ggplot(Imp_Q_micro_filt2, aes(x=Overall, y=Taxonomy)) + geom_point(stat = "identity") + 
geom_vline(xintercept=50, color="tomato") + theme(axis.text=element_text(color="black")) + theme_bw()  + 
scale_x_reverse() +ggtitle("Microbial feature ranking, paired biopsy classifier") + xlab("Overall feature importance")


df[1:2,]



