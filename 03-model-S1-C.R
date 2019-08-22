rm(list=ls())
load('./data/raw.data.RData')
y.ga <- raw.data$train$phenotype$GA
names(y.ga) <- raw.data$train$phenotype$SampleID
x.ga <- raw.data$train$x[,names(y.ga)]
test.data <- t(raw.data$test$x)
rm(raw.data)
library('glmnet')
x.train <- t(x.ga)
y.train <- log2(y.ga)
set.seed(768)
fold10 <- caret::createFolds(
  y=y.train,
  k=10)
fold.model <- list()
rmse <- c()
for(i in 1:length(fold10)){
  cat(paste0('Cross-validation: i=',i,'\n'))
  train.set <- names(y.train[unlist(fold10[-i])])
  test.set <- names(y.train[unlist(fold10[i])])
  tr.x <- x.train[train.set,]
  tr.y <- y.train[train.set]
  a1.fit <- cv.glmnet(
    tr.x,
    tr.y,
    type.measure = 'mse',
    alpha=1,
    family='gaussian'
  )
  tt.x <- x.train[test.set,]
  tt.y <- y.train[test.set]
  a1.pred <- predict(
    a1.fit,
    s=a1.fit$lambda.1se,
    new=tt.x
  )
  rmse[i] <- RMSE(a1.pred,tt.y)
  fold.model[[i]] <- a1.fit
}
lasso.mods <- fold.model
save(rmse,lasso.mods,test.data,
     file='./data/models/model_S1_C.RData')
#----------------------
# PREDICTION
#----------------------
rm(list=ls())
load('./data/models/model_S1_C.RData')
pdf('./figures/03-Model-S1-C-performance-median.pdf',
    width = 4,height = 5)
hist(2^rmse,
     col = '#CBC9E2',
     xlab = 'RMSE',main='LASSO | 10 fold CV')
dev.off()
pred.x <- list()
for(i in 1:length(lasso.mods)){
  a1.fit <- lasso.mods[[i]]
  pred.x[[i]] <- predict(
    a1.fit,
    s=a1.fit$lambda.1se,
    new=test.data
  )
}
xx <- do.call('cbind',pred.x)
pred.test.x <-2^apply(xx,1,mean)
out <- data.frame(SampleID = names(pred.test.x),stringsAsFactors = FALSE)
out$GA <- pred.test.x
rownames(out) <- out$SampleID
test.sample <- read.csv(
  './data/TeamX_SC1_prediction_MeanModel.csv',
  stringsAsFactors = FALSE
)
write.csv(
  out[test.sample$SampleID,],
  './data/prediction/TeamIGIB_SC1_MODS1C.csv',
  row.names = FALSE,
  quote = FALSE
)
