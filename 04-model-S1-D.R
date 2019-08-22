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
loocv.bias.rse <- list()
loocv.bias.test.pred <- list()
y.sort <- sort(y.train)
# 3% top and bottom
n3TB <- round((length(y.sort)/100)*3)
loocv.bias.loop <- setdiff(
  names(y.sort),
  c(head(names(y.sort),n3TB),tail(names(y.sort),n3TB))
)
dir.create('./data/models/loocv-bias/',showWarnings = FALSE)
loocv.bias.rmse <- list()
loocv.bias.test.pred <- list()
for(i in 1:length(loocv.bias.loop)){
  print(i)
  train.set <- setdiff(names(y.train),loocv.bias.loop[i])
  test.set <- loocv.bias.loop[i]
  tr.x <- x.train[train.set,]
  tr.y <- y.train[train.set]
  a1.fit <- cv.glmnet(
    tr.x,
    tr.y,
    type.measure = 'mse',
    alpha=1,
    family='gaussian'
  )
  a0.fit <- cv.glmnet(
    tr.x,
    tr.y,
    type.measure = 'mse',
    alpha=0,
    family='gaussian'
  )
  tt.x <- t(as.matrix(x.train[test.set,]))
  tt.y <- y.train[test.set]
  a1.pred <- predict(
    a1.fit,
    s=a1.fit$lambda.1se,
    new=tt.x
  )
  a0.pred <- predict(
    a0.fit,
    s=a0.fit$lambda.1se,
    new=tt.x
  )
  loocv.bias.rmse[[i]] <- c(
    'a1' = caret::RMSE(tt.y, a1.pred),
    'a0' = caret::RMSE(tt.y, a0.pred)
  )
  loocv.bias.model <- list(
    'a1' = a1.fit,
    'a0' = a0.fit
  )
  loocv.bias.test.pred[[i]] <- list(
    'pred.a1' = predict(
      a1.fit,
      s=a1.fit$lambda.1se,
      new=test.data
    ),
    'pred.a0' = predict(
      a0.fit,
      s=a0.fit$lambda.1se,
      new=test.data
    )
  )
  mod.name <- paste0(
    './data/models/loocv-bias/loocv-',
    test.set,
    sprintf("%03d", i),
    '.RData'
  )
  save(
    loocv.bias.model,
    file = mod.name
  )
}
save(loocv.bias.test.pred,
     loocv.bias.rmse,
     file = './data/models/loocv.bias.outputs.RData'
)
#-----------------------
rm(list=ls())
load('./data/models/loocv.bias.outputs.RData')

a1.out <- lapply(loocv.bias.test.pred,
                 function(x){
                   x$pred.a1
                 })
xx <- do.call('cbind',a1.out)
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
  './data/prediction/TeamIGIB_SC1_MODS1D_loocv_bias.csv',
  row.names = FALSE,
  quote = FALSE
)
