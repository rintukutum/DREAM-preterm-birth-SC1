#-----------------------
dir.create('./data/',showWarnings = FALSE)
dir.create('./figures/',showWarnings = FALSE)
dir.create('./data/prediction/',showWarnings = FALSE)
dir.create('./data/models/',showWarnings = FALSE)
#-----------------------
# Download files
# 1. TeamX_SC1_prediction_MeanModel.csv
# 2. anoSC1_v11_nokey.csv
# 3. HTA20_RMA.RData
#-----------------------
rm(list=ls())
# write a code to download the data from Synapse
pheno <- read.csv(
  './data/anoSC1_v11_nokey.csv',
  stringsAsFactors = FALSE
)
pheno.train <- pheno[pheno$Train == 1,]
pheno.test <- pheno[pheno$Train == 0,]

load('./data/HTA20_RMA.RData')

x.train <- eset_HTA20[,pheno.train$SampleID]
x.test <- eset_HTA20[,pheno.test$SampleID]
# remove from memory
rm(eset_HTA20)

raw.data <- list(
  train = list(
    phenotype = pheno.train,
    x = x.train
  ),
  test = list(
    phenotype = pheno.test,
    x = x.test
  )
)
save(raw.data, file = './data/raw.data.RData')