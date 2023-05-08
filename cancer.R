# Load packages
library(corrplot)
library(VIM)
library(ggplot2)
library(gridExtra)
library(scales)
library(tidyr)
library(dplyr)
library(readr)
library(ggformula)
library(caret)
library(rpart) 
library(rpart.plot)
library(kernlab)
library(randomForest)
library(car)
library(glmnet)
library(xgboost)


# Load data
df <- read.csv("./METABRIC_RNA_Mutation.csv", header=TRUE, sep=",", stringsAsFactors = FALSE);
dim(df)

## Preparing and Cleaning the Data

# To replace “ ” with  NA
df[df == ""] <- NA            
df_t = df[, c(1:31)] 

# To calculate percent of missing data for each variable
M_Val <- function(x){mm_v=round((sum(is.na(x))/length(x))*100,3)} 
nam_1 <- colnames(df_t)
c <- apply(df_t,2,M_Val)

df_p <- data.frame(nam_1,c)
rownames(df_p)<-NULL
df_p1 <- df_p %>%  filter(c > 5 )  %>%  mutate(c1=c-2)
gf_boxplot(c~nam_1 , data = df_p) %>%   
  gf_labs(title = "Missing data percentage",
          x = "",  
          y = "Percent (%)")%>%   gf_theme(axis.text.x=element_text(angle=65, hjust=1)) %>%
  gf_label(c1~nam_1, 
           label =~ c, 
           data = df_p1, col = "black", 
           fill =~ nam_1, alpha = .2,
           show.legend = FALSE)


# To remove Redundant  column
df <- df %>% 
  select(-patient_id,-cancer_type,-death_from_cancer,-cancer_type_detailed,
         -oncotree_code,-cohort,-tumor_stage,-overall_survival_months,-primary_tumor_laterality,
         -X3.gene_classifier_subtype,-er_status_measured_by_ihc,-her2_status_measured_by_snp6)

clinical = df[, c(1:19)]    # To make clinical data set

# To remove factor with one attribute
clinical <- clinical[!( clinical$tumor_other_histologic_subtype=="Metaplastic"| clinical$pam50_._claudin.low_subtype=="NC"),]
clinical$tumor_other_histologic_subtype <- as.factor(clinical$tumor_other_histologic_subtype); 
clinical$pam50_._claudin.low_subtype <- as.factor(clinical$pam50_._claudin.low_subtype); 

# To change character variables to factor 
clinical$type_of_breast_surgery <- as.factor(clinical$type_of_breast_surgery)
clinical$cellularity <- as.factor(clinical$cellularity); 
clinical$chemotherapy <- as.factor(clinical$chemotherapy)
clinical$pam50_._claudin.low_subtype <- as.factor(clinical$pam50_._claudin.low_subtype); 
clinical$er_status <- as.factor(clinical$er_status); 
clinical$neoplasm_histologic_grade <- as.factor(clinical$neoplasm_histologic_grade);
clinical$her2_status <- as.factor(clinical$her2_status);
clinical$tumor_other_histologic_subtype <- as.factor(clinical$tumor_other_histologic_subtype)
clinical$hormone_therapy <- as.factor(clinical$hormone_therapy);
clinical$inferred_menopausal_state<- as.factor(clinical$inferred_menopausal_state)
clinical$integrative_cluster<- as.factor(clinical$integrative_cluster ); 
#clinical$primary_tumor_laterality<- as.factor(clinical$primary_tumor_laterality)
clinical$pr_status<- as.factor(clinical$pr_status); 
clinical$radio_therapy <- as.factor(clinical$radio_therapy)
clinical <- clinical %>%  mutate(lymph_positive = ifelse(lymph_nodes_examined_positive == 0, "0", "1"))
clinical$lymph_positive <- as.factor(clinical$lymph_positive); clinical$lymph_nodes_examined_positive <- NULL
clinical$overall_survival<- as.factor(clinical$overall_survival); 

# Drop the missing values
clinical <- na.omit(clinical)


# Correlation plot
clinical_n <- select_if(clinical, is.numeric) 
corrmatrix <- cor(clinical_n)
corrplot::corrplot(corrmatrix, method="shade", type = "upper",tl.cex=.6	
                   , tl.col="black", title="Correlation Plot",number.font = 2, mar=c(0,0,1,0))

# Exploring numerical data
par(mfrow=c(4 ,2))
clinical %>% gf_boxplot(age_at_diagnosis~overall_survival) %>%   
  gf_labs(x = "Survival",  
          y = "Age") %>%   
  gf_theme(axis.text = element_text(size = 12),
           axis.title = element_text(size = 12),
           legend.title=element_text(size=14), 
           legend.text=element_text(size=9))
QQlabels = c("0" = "overall_survival=0","1" = "overall_survival=1")
ggplot(data = clinical, aes(sample=age_at_diagnosis)) + 
  stat_qq() + 
  stat_qqline() + 
  facet_grid(.~overall_survival,labeller=labeller(overall_survival=QQlabels))

clinical %>% gf_boxplot(nottingham_prognostic_index~overall_survival) %>%   
  gf_labs(x = "Survival",  
          y = "nottingham_prognostic_index") %>%   
  gf_theme(axis.text = element_text(size = 12),
           axis.title = element_text(size = 12),
           legend.title=element_text(size=14), 
           legend.text=element_text(size=9))
QQlabels = c("0" = "overall_survival=0","1" = "overall_survival=1")
ggplot(data = clinical, aes(sample=nottingham_prognostic_index)) + 
  stat_qq() + 
  stat_qqline() + 
  facet_grid(.~overall_survival,labeller=labeller(overall_survival=QQlabels))

clinical %>% gf_boxplot(mutation_count~overall_survival) %>%   
  gf_labs(x = "Survival",  
          y = "Mutation") %>%   
  gf_theme(axis.text = element_text(size = 12),
           axis.title = element_text(size = 12),
           legend.title=element_text(size=14), 
           legend.text=element_text(size=9))
QQlabels = c("0" = "overall_survival=0","1" = "overall_survival=1")
ggplot(data = clinical, aes(sample=mutation_count)) + 
  stat_qq() + 
  stat_qqline() + 
  facet_grid(.~overall_survival,labeller=labeller(overall_survival=QQlabels))

clinical %>% gf_boxplot(tumor_size~overall_survival) %>%   
  gf_labs(x = "Survival",  
          y = "Tumor size") %>%   
  gf_theme(axis.text = element_text(size = 12),
           axis.title = element_text(size = 12),
           legend.title=element_text(size=14), 
           legend.text=element_text(size=9))
QQlabels = c("0" = "overall_survival=0","1" = "overall_survival=1")
ggplot(data = clinical, aes(sample=tumor_size)) + 
  stat_qq() + 
  stat_qqline() + 
  facet_grid(.~overall_survival,labeller=labeller(overall_survival=QQlabels))

# The logarithm transformation
clinical$tumor_size_T <- log10((clinical$tumor_size)); clinical$tumor_size <- NULL
clinical$mutation_count_T <- log10((clinical$mutation_count)); clinical$mutation_count <- NULL

# Plot data
par(mfrow=c(2 ,2))
hist(clinical_n$tumor_size , main = "", xlab = "Tumor size", freq = FALSE )
curve(dnorm(x, mean = mean(clinical_n$tumor_size), sd = sd(clinical_n$tumor_size)),col = "red", add = TRUE)

hist(clinical$tumor_size_T , main = "", xlab = "Tumor size ", freq = FALSE )
curve(dnorm(x, mean = mean(clinical$tumor_size_T), sd = sd(clinical$tumor_size_T)),col = "red", add = TRUE)

hist(clinical_n$mutation_count, main = "", xlab = "Mutation count", freq = FALSE )
curve(dnorm(x, mean = mean(clinical_n$mutation_count), sd = sd(clinical_n$mutation_count)),col = "red", add = TRUE)

hist(clinical$mutation_count_T, main = "", xlab = "Mutation count", freq = FALSE )
curve(dnorm(x, mean = mean(clinical$mutation_count_T), sd = sd(clinical$mutation_count_T)),col = "red", add = TRUE)


# Rename levels in response   
clinical_e <- clinical
make.names("overall_survival", unique = TRUE, allow_ = FALSE); 
clinical_e  %>% 
  mutate(overall_survival = factor(overall_survival, 
                                   labels = make.names(levels(overall_survival))));
levels(clinical_e$overall_survival) <- c("D", "L")


## Training models
# Spliting samples
set.seed(100)
train_indices <- sample(nrow(clinical_e),nrow(clinical_e)*0.9,replace = F)
train <- clinical_e[train_indices,]
test <- clinical_e[-train_indices,]
y <- test$overall_survival

# ROC curve
par(mfrow=c(1 ,1))
ROC_AUC <- function(pred, gold){
  gold <- factor(gold, levels = c("0","1"))
  TPR_list <- c(0)
  FPR_list <- c(0)
  AUC <- 0
  for (theta in seq(1,0,-0.01)) {
    pred_c <- factor(1*(pred>theta), levels = c("0","1"))
    cm <- table(gold,pred_c)
    TPR <- cm[2,2]/(cm[2,1]+cm[2,2])
    FPR <- cm[1,2]/(cm[1,1]+cm[1,2])
    AUC <- AUC + (FPR-FPR_list[length(FPR_list)])*(TPR+TPR_list[length(TPR_list)])/2
    TPR_list <- c(TPR_list,TPR)
    FPR_list <- c(FPR_list,FPR)
  }
  return(list(AUC=AUC, TPR=TPR_list[-1], FPR=FPR_list[-1]))
}

# 1. logistic
logistic_model <- glm(I(overall_survival=='L')~., data=train, family = 'binomial')
summary(logistic_model) # regression coefficients
probs <- predict.glm(logistic_model, newdata = test, type = 'response')
pred <- 1*(probs>0.5)
pred <- factor(pred, labels = c('D','L'))
table(y, pred)
mean(y==pred) # acc
r1 <- ROC_AUC(probs, 1*(y=='L'))


# 2. lasso
train_X <- model.matrix(~.-overall_survival, data = train)
train_y <- 1*(train$overall_survival=='L')
test_X <- model.matrix(~.-overall_survival, data = test)

cv_lasso <- cv.glmnet(x = train_X, y = train_y,type.measure = 'class',nfolds = 10,family='binomial',gamma = 1)
lasso_model <- glmnet(x = train_X, y = train_y, family = 'binomial', alpha = 1, lambda = cv_lasso$lambda.min)
lasso_model$beta

probs <- predict(lasso_model, newx = test_X, type = 'response')
pred <- factor(1*(probs>0.5), labels = c('D','L'))
table(y, pred)
mean(y==pred) # acc
r2 <- ROC_AUC(probs, 1*(y=='L'))


# 3. ridge
cv_ridge <- cv.glmnet(x = train_X, y = train_y,type.measure = 'class',nfolds = 10,family='binomial',gamma = 0)
ridge_model <- glmnet(x = train_X, y = train_y, family = 'binomial', alpha = 0, lambda = cv_lasso$lambda.min)
ridge_model$beta

probs <- predict(ridge_model, newx = test_X, type = 'response')
pred <- factor(1*(probs>0.5), labels = c('D','L'))
table(y, pred)
mean(y==pred) # acc
r3 <- ROC_AUC(probs, 1*(y=='L'))


# 4. elastic net
cv.elnet <- cv.glmnet(x = train_X, y = train_y, type.measure = 'class',nfolds = 10,family='binomial')
coef(cv.elnet)
probs <- predict(cv.elnet, newx=test_X, s='lambda.min', gamma='gamma.min', type='response')
pred <- factor(1*(probs>0.5), labels = c('D','L'))
table(y, pred)
mean(y==pred) # acc
r4 <- ROC_AUC(probs, 1*(y=='L'))


# 5. decision tree CART
dtree <- rpart(overall_survival~., data=train, method = 'class', parms=list(split="gini"))
tree <- prune(dtree,cp=dtree$cptable[which.min(dtree$cptable[,"xerror"]),"CP"])
tree$variable.importance
# par(mfrow=c(1,1))
rpart.plot(tree,branch=1, type=4,fallen.leaves=T,cex=0.8)
pred <- predict(tree, newdata=test, type='class')
probs <- predict(tree, newdata=test, type='prob')[,2]
table(y, pred)
mean(y==pred) # acc
r5 <- ROC_AUC(probs, 1*(y=='L'))


# 6. decision tree ID3
dtree <- rpart(overall_survival~., data=train, method = 'class', parms=list(split="information"))
tree <- prune(dtree,cp=dtree$cptable[which.min(dtree$cptable[,"xerror"]),"CP"])
tree$variable.importance
rpart.plot(tree,branch=1, type=4,fallen.leaves=T,cex=0.8)
pred <- predict(tree, newdata=test, type='class')
probs <- predict(tree, newdata=test, type='prob')[,2]
table(y, pred)
mean(y==pred) # acc
r6 <- ROC_AUC(probs, 1*(y=='L'))


# 7. random forest
rf <- randomForest(overall_survival~., data=train, importance=TRUE)
print(rf)
rf$importance # importance of variables
pred <- predict(rf, test)
probs <- predict(rf, newdata=test, type='prob')[,2]
table(y, pred)
mean(y==pred) # acc
r7 <- ROC_AUC(probs, 1*(y=='L'))


# 8. xgboost
library(Matrix)
train_matrix <- sparse.model.matrix(overall_survival ~ .-1, data = train)
test_matrix <- sparse.model.matrix(overall_survival ~ .-1, data = test)
train_y <- as.numeric(train$overall_survival=='L')
test_y <-  as.numeric(test$overall_survival=='L')
dtrain <- xgb.DMatrix(data=train_matrix,label=train_y)
dtest <- xgb.DMatrix(data=test_matrix,label=test_y)
xgb <- xgboost(data=dtrain,max_depth=6, eta=0.5,  objective='binary:logistic', nround=25) 
# importance of variables
importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)  
print(importance)
probs <- predict(xgb,newdata = dtest)
pred <- factor(1*(probs>0.5), labels = c('D','L'))
table(y, pred)
mean(y==pred) # acc
r8 <- ROC_AUC(probs, 1*(y=='L'))

# ROC plot
plot(r1$FPR,r1$TPR,type = 'l',col = 'blue',xlab = "FPR",ylab = "TPR",legend = r1$AUC)
lines(r2$FPR,r2$TPR,col = 'red')
lines(r3$FPR,r3$TPR,col = 'green')
lines(r4$FPR,r4$TPR,col = 'purple')
lines(r5$FPR,r5$TPR,col = 'pink')
lines(r6$FPR,r6$TPR,col = 'yellow')
lines(r7$FPR,r7$TPR,col = 'brown')
lines(r8$FPR,r8$TPR,col = 'gray')
abline(a=0, b=1, lty=2)
legend('bottomright',legend = c(paste0('logistic:',r1$AUC),paste0('lasso:',r2$AUC),
       paste0('ridge:',r3$AUC), paste0('elastic net:',r4$AUC), paste0('tree CART:',r5$AUC),
        paste0('tree ID3:',r6$AUC), paste0('random forest:',r7$AUC), paste0('xgboost:',r8$AUC)),
       col = c('blue','red','green','purple','pink','yellow','brown','gray'), lty = 1,cex = 0.5)
title('ROC curve')





