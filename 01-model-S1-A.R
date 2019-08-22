rm(list=ls())
#----------------------------
# PREPARE FOR MODEL
#----------------------------
load('./data/raw.data.RData')
y.ga <- raw.data$train$phenotype$GA
names(y.ga) <- raw.data$train$phenotype$SampleID
x.ga <- raw.data$train$x[,names(y.ga)]
test.data <- raw.data$test
rm(raw.data)
#----------------------------
# FEATURE REDUCTION
# 2-5 mins
#----------------------------
library('sparsepca')
out <- spca(x.ga, k=3, alpha=1e-3, beta=1e-3, center = TRUE, scale = FALSE, verbose=0)
idx.pc1 <- which(out$loadings[,1] != 0)
idx.pc2 <- which(out$loadings[,2] != 0)
idx.pc3 <- which(out$loadings[,3] != 0)
feature.idx <- c(idx.pc1,idx.pc2,idx.pc3)
x.ga.model <- x.ga[feature.idx,]
rm(x.ga)
#----------------------------
# TRAIN DATA
#----------------------------
x.train <- t(x.ga.model)
rm(x.ga.model)
y.train <- y.ga
rm(y.ga)
featureNames <- colnames(x.train)
names(featureNames) <- paste0('F',1:length(featureNames))
colnames(x.train) <- names(featureNames)
test.xx <- test.data$x[as.character(featureNames),]
rownames(test.xx) <- names(featureNames)
test.xxx <- t(test.xx)
x.train <- as.matrix(x.train)
#----------------------------
# MODEL
#----------------------------

#----------------------------
# FEATURE SELECTION WITH Boruta
# 5-10 mins
#----------------------------
library(Boruta)
set.seed(7853)
bt <- Boruta(
  x = x.train,
  y=y.train,
  maxRuns = 500
)
idx.bt.var <- bt$finalDecision %in% c('Confirmed','Tentative')
x.vars <- names(bt$finalDecision[idx.bt.var])
#----------------------------
# 10 fold cross-validation
#----------------------------
library(caret)
set.seed(768)
fold10 <- createFolds(y=y.train,k=10)
#----------------------------
# svm and random forest
#----------------------------
library(randomForest)
library(e1071)
fold.model <- list()
cor <- list()
for(i in 1:length(fold10)){
  cat(paste0('Cross-validation: i=',i,'\n'))
  train.set <- names(y.train[unlist(fold10[-i])])
  test.set <- names(y.train[unlist(fold10[i])])
  tr.x <- x.train[train.set,x.vars]
  tr.y <- y.train[train.set]
  tt.x <- x.train[test.set,x.vars]
  tt.y <- y.train[test.set]
  set.seed(7686+i)
  rf <- randomForest(
    x = tr.x,
    y = tr.y,
    ntree = 2000
  )
  pred.rf <- predict(rf,tt.x)
  rf.cor <- cor(pred.rf,tt.y)
  sv.rad <- svm(x=tr.x,y=tr.y,kernel = 'radial')
  sv.li <- svm(x=tr.x,y=tr.y,kernel = 'linear')
  sv.pol <- svm(x=tr.x,y=tr.y,kernel = 'polynomial')
  sv.sig <- svm(x=tr.x,y=tr.y,kernel = 'sigmoid')
  pred.sv <- list(
    'linear' = predict(sv.li,tt.x),
    'polynomial' = predict(sv.pol,tt.x),
    'radial' = predict(sv.rad,tt.x),
    'sigmoid' = predict(sv.sig,tt.x)
  )
  sv.cor <- sapply(
    pred.sv,
    function(x){
      cor(
        x,
        tt.y)})
  cor[[i]] <- c(svm=sv.cor,rf=rf.cor)
  fold.model[[i]] <- list(
    rf=rf,
    svm=list(
      'linear' = sv.li,
      'polynomial' = sv.pol,
      'radial' = sv.rad,
      'sigmoid' = sv.sig
    ))
}
save(fold.model,cor,test.xxx,x.vars,y.train,
     file = './data/models/model_S1_A.RData')
#----------------------------
# TESTING
#----------------------------
rm(list=ls())
load('./data/models/model_S1_A.RData')
rf.test <- list()
for(i in 1:length(fold.model)){
  rf.test[[i]] <- predict(fold.model[[i]]$rf,test.xxx[,x.vars])
}

svm.linear.test <- list()
for(i in 1:length(fold.model)){
  svm.linear.test[[i]] <- predict(
    fold.model[[i]]$svm$linear,
    test.xxx[,x.vars])
}
svm.poly.test <- list()
for(i in 1:length(fold.model)){
  svm.poly.test[[i]] <- predict(
    fold.model[[i]]$svm$polynomial,
    test.xxx[,x.vars])
}

svm.rad.test <- list()
for(i in 1:length(fold.model)){
  svm.rad.test[[i]] <- predict(
    fold.model[[i]]$svm$radial,
    test.xxx[,x.vars])
}

svm.sig.test <- list()
for(i in 1:length(fold.model)){
  svm.sig.test[[i]] <- predict(
    fold.model[[i]]$svm$sigmoid,
    test.xxx[,x.vars])
}
svm.sig.df <- do.call('cbind',svm.sig.test)
svm.rad.df <- do.call('cbind',svm.rad.test)
svm.lin.df <- do.call('cbind',svm.linear.test)
svm.pol.df <- do.call('cbind',svm.poly.test)
rf.test.df <- do.call('cbind',rf.test)
#----------------------------
# FIGURE
#----------------------------
pdf('./figures/01-Model-S1-A-performance-mdeian.pdf',
    width = 4,height = 5)
boxplot(list(
  Original = y.train,
  SVM.Rad = apply(svm.rad.df,1,median),
  SVM.Sig = apply(svm.sig.df,1,median),
  SVM.Lin = apply(svm.lin.df,1,median),
  SVM.Pol = apply(svm.pol.df,1,median),
  RF = apply(rf.test.df,1,median)
  ),
  col = c(
    Original = "#4DAF4A",
    SVM.Rad = "#EFF3FF",
    SVM.Sig = "#BDD7E7",
    SVM.Lin = "#6BAED6",
    SVM.Pol = "#2171B5",
    RF = "#CBC9E2"),
  las =2,
  ylab = 'Gestational Age(weeks)',
  main ='Model-S1-A performance (10 fold CV)'
)
dev.off()
#----------------------------
# SVM radial median
#----------------------------
svm.rad.out <- data.frame(
  apply(svm.rad.df,1,median)
)
colnames(svm.rad.out)[1] <- 'GA'
svm.rad.out$SampleID <- rownames(svm.rad.out)
test.sample <- read.csv(
  './data/TeamX_SC1_prediction_MeanModel.csv',
  stringsAsFactors = FALSE
)
dir.create('./data/prediction/',showWarnings = FALSE)

#----------------------------
#  PREDICTION RESULT
#----------------------------
write.csv(
  svm.rad.out[test.sample$SampleID,2:1],
  './data/prediction/TeamIGIB_SC1_MODS1A.csv',
  row.names = FALSE,
  quote = FALSE
)