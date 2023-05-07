set.seed(5000)

library(rpart)
library(rpart.plot)
library(AUC)
library(pROC)
library(party)
library(partykit)
library(ROCR)
library(evtree)
library(tree)
library(caret)
library(C50)
library(ggplot2)
library(BiocManager)
library(svMisc)
library(glm2)

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))


dataset <- readRDS("datasets/GSE102299-smartseq2.rds")
dataset <- updateObject(dataset)
colData(dataset)[,14]
head(colData(dataset))
table(colData(dataset)[,14])

(dataset <- experiments(dataset)[["gene"]])
dataset <- assays(dataset)[["avetxlength"]]
dataset <- t(dataset)
dim(dataset)

#Based on characteristics_ch1.4 factor (treatment: IL25 & treatment: IL33) #
dataset <- dataset[1:470,]
dim(dataset)
dataset[1:3,1:10]
row.names(dataset)[1:188] <- "1"
row.names(dataset)[189:470] <- "0"
dataset[1:3,1:10]

write.csv(dataset, "Test.csv")
dataset <- read.csv("Test.csv")
dataset<- data.frame(dataset)
#-------------------------------------------------------------------------
gene_id <- c()
p_values <- c()
for(n in 2:ncol(dataset)) {
  wilcox <- wilcox.test(dataset[[n]] ~ dataset[[1]], data = dataset, exact = T, paired = FALSE)
  p_values[n-1] <- wilcox$p.value
  cat("\n", n)
}

gene_id <- c(1:(ncol(dataset)-1))
head(gene_id)
head(p_values)
lowest_p <- data.frame(gene_id, p_values)
lowest_p <- lowest_p[order(p_values),]
head(lowest_p)
tail(lowest_p)
candidate_genes <- lowest_p$gene_id
candidate_genes <- candidate_genes[1:1000]
head(candidate_genes)
dataset <- dataset[ , c(1, candidate_genes)]

glm_runtime <- c()
glm_precision <- c()
glm_recall <- c()
glm_F1score <- c()
glm_auc <- c()

for(k in 1:100) {
  cat("\014")
  cat(k)
  cat("\n")
  
  shuffled <- dataset[sample(nrow(dataset), replace=FALSE), ] # Shuffling dataset
  dim(dataset)
  dim(shuffled)
  
  glm_runtime_avearge <- c()
  glm_precision_avearge <- c()
  glm_recall_avearge <- c()
  glm_F1score_avearge <- c()
  glm_auc_avearge <- c()

  
  i=1 
  j=47
  for(x in 1:10) {
    test_set <- shuffled[i:j, ]
    train_set <- shuffled[-c(i:j), ]
    #--------------------------------------------------------------------------

    start_time <- Sys.time()
    glm_model <- glm(X ~., data = train_set, family = "binomial")
    end_time <- Sys.time()
    glm_runtime_avearge <- append(glm_runtime_avearge, (as.numeric(difftime(end_time,start_time, units = c("secs")))))
    result <- predict(glm_model, test_set)
    
    confu_matrix <- table(factor(test_set$X, levels = c(0,1)), factor(result >= 0.5,levels = c(F,T)))
    
    stopifnot((ncol(confu_matrix))&&(nrow(confu_matrix)) == 2)
    
    glm_precision_avearge <- append(glm_precision_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))
    glm_recall_avearge <- append(glm_recall_avearge, confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))
    glm_F1score_avearge <- append(glm_F1score_avearge, (2*((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))*(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1])))/((confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[1,2]))+(confu_matrix[2,2]/(confu_matrix[2,2]+confu_matrix[2,1]))))
    )
    
    test_error <- try({roc <- roc(result, test_set$X)}, silent = T)
    if(class(test_error)=="try-error") {
      glm_auc_avearge <- append(glm_auc_avearge, 1)
    } else {
      roc <- roc(result, test_set$X)
      glm_auc_avearge <- append(glm_auc_avearge, roc$auc[1])
    }
    test_error <- NULL
    #--------------------------------------------------------------------------
    cat("\014")
    cat(k)
    cat("\n")
 
    
    i=j+1
    if(i == 424) {
      j=j+47
    } else {
      j=j+47
    }
  }

  
  glm_runtime <- append(glm_runtime, mean(glm_runtime_avearge))
  glm_precision <- append(glm_precision, mean(glm_precision_avearge))
  glm_recall <- append(glm_recall, mean(glm_recall_avearge))
  glm_F1score <- append(glm_F1score, mean(glm_F1score_avearge))
  glm_auc <- append(glm_auc, mean(glm_auc_avearge))
  
}

Runtime_table <- data.frame("Runtime" = glm_runtime)
write.csv(Runtime_table, "~/Desktop/glm()_project/GSE102299-smartseq2/Runtime_table.csv")
pdf("~/Desktop/glm()_project/GSE102299-smartseq2/Runtime.pdf")
boxplot(Runtime_table$Runtime, outline=FALSE, main = "Runtime", xlab = "GSE102299-smartseq2", ylab = "Seconds", col=c("#800020")) #outline=FALSE to remove extreme values
text(y=fivenum(Runtime_table$Runtime), labels = round(fivenum(Runtime_table$Runtime), digits = 3), x=0.74)
legend("topright", bty = "n", legend = summary(Runtime_table))
dev.off()

Precision_table <- data.frame("Precision" = glm_precision)
write.csv(Precision_table, "~/Desktop/glm()_project/GSE102299-smartseq2/Precision_table.csv")
pdf("~/Desktop/glm()_project/GSE102299-smartseq2/Precision.pdf")
boxplot(Precision_table$Precision, outline=FALSE, main = "Precision", xlab = "GSE102299-smartseq2", ylab = "Score", col=c("#800020")) #outline=FALSE to remove extreme values
text(y=fivenum(Precision_table$Precision), labels = round(fivenum(Precision_table$Precision), digits = 3), x=0.74)
legend("topright", bty = "n", legend = summary(Precision_table))
dev.off()

Recall_table <- data.frame("Recall" = glm_recall)
write.csv(Recall_table, "~/Desktop/glm()_project/GSE102299-smartseq2/Recall_table.csv")
pdf("~/Desktop/glm()_project/GSE102299-smartseq2/Recall.pdf")
boxplot(Recall_table$Recall, outline=FALSE, main = "Recall", xlab = "GSE102299-smartseq2", ylab = "Score", col=c("#800020")) #outline=FALSE to remove extreme values
text(y=fivenum(Recall_table$Recall), labels = round(fivenum(Recall_table$Recall), digits = 3), x=0.74)
legend("topright", bty = "n", legend = summary(Recall_table))
dev.off()

F1score_table <- data.frame("F1score" = glm_F1score)
write.csv(F1score_table, "~/Desktop/glm()_project/GSE102299-smartseq2/F1score_table.csv")
pdf("~/Desktop/glm()_project/GSE102299-smartseq2/F1score.pdf")
boxplot(F1score_table$F1score, outline=FALSE, main = "F1score", xlab = "GSE102299-smartseq2", ylab = "Score", col=c("#800020")) #outline=FALSE to remove extreme values
text(y=fivenum(F1score_table$F1score), labels = round(fivenum(F1score_table$F1score), digits = 3), x=0.74)
legend("topright", bty = "n", legend = summary(F1score_table))
dev.off()

AUC_table <- data.frame("AUC" = glm_auc)
write.csv(AUC_table, "~/Desktop/glm()_project/GSE102299-smartseq2/AUC_table.csv")
pdf("~/Desktop/glm()_project/GSE102299-smartseq2/AUC.pdf")
boxplot(AUC_table$AUC, outline=FALSE, main = "AUC", xlab = "GSE102299-smartseq2", ylab = "Score", col=c("#800020")) #outline=FALSE to remove extreme values
text(y=fivenum(AUC_table$AUC), labels = round(fivenum(AUC_table$AUC), digits = 3), x=0.74)
legend("topright", bty = "n", legend = summary(AUC_table))
dev.off()
