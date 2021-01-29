require('MLSeq')

# --- NEED RAW COUNT DATA  --- #

filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- read.table(filepath, header=TRUE)
class <-S4Vectors::DataFrame(condition = factor(rep(c("N","T"), c(29, 29))))

# load in raw data
transcripts_to_genes <- read_csv("data/human_xy_rnaseq.csv")

loadToy <- function(idx){
  toy1 <- read_csv(sprintf("data/toy_counts/toy_df_%s.csv", idx)) # --> filter for XY chromosome only
  toy1xy <- toy1 %>% 
    semi_join(transcripts_to_genes, by=c("gene_name"="transcript")) %>% 
    data.frame()
  dim(toy1xy)
  rownames(toy1xy) <- toy1xy$gene_name
  toy1xy$gene_name <- NULL
  return(toy1xy)
}

toy_df <- do.call(cbind, lapply(1:5, loadToy))
toy_df2 <- toy_df[,unique(colnames(toy_df))]

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdld") # get some class labels
keep_cols <- comb_metadata %>% filter(sample_acc %in% colnames(toy_df2))
sex_lab <- keep_cols %>% pull(metadata_sex)
names(sex_lab) <- keep_cols %>% pull(sample_acc)
sex_lab2 <- sex_lab[colnames(toy_df2)]
table(sex_lab2)
class <-S4Vectors::DataFrame(condition = factor(sex_lab2))

# put together a reasonably sized df
# try playing around with the transcript-level data

  
library(DESeq2)
set.seed(2128)
# select top 100 features
vars <- sort(apply(toy_df2, 1, var, na.rm = TRUE), decreasing = TRUE)
data <- apply(toy_df2[names(vars)[1:1000], ], c(1,2), round)

nTest <- ceiling(ncol(data) * 0.3)
ind <- sample(ncol(data), nTest, FALSE)
# Minimum count is set to 1 in order to prevent 0 division problem within
# classification models.
data.train <- as.matrix(data[ ,-ind] + 1)
data.test <- as.matrix(data[ ,ind] + 1)
classtr <- DataFrame(condition = class[-ind, ])
classts <- DataFrame(condition = class[ind, ])

data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = classtr,
                                      design = formula(~condition))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = classts,
                                     design = formula(~condition))

# ---- normalization ---- #
#  variations: library-size, majority fragments?, within-sample (gene length, GC content)
#
# "deseq median ratio normalization"
#   - divide each sample by the mean of transcript counts
# "trimmed mean of M values" (TMM)
#   - trims the data in upper log-fold changes to minimize btw-sample / intensity changes

# can use discrete classifiers after this (PLDA, NBLDA)
# -OR-
# transform --> similar to microarrays
#  - log-cpm (less skewed but unequal genewise variances)
#  - control w/ vst (variance stabilizing transformation), rlog (regularized log transform)
#      -rlog similar, but shrunken for genes with lower counts

#alpha=seq(0,1, 0.1)
#lambda=c(0.01, 0.02, 0.05,  0.10, 0.20, 0.50)
#tune_df <- cross_df(list("alpha"=alpha, "lambda"=lambda))


fit <- classify(data=data.trainS4, method="glmnet",
                preProcessing="deseq-vst", ref="male",
                tuneLength=10,
                control = trainControl(method = "repeatedcv", 
                                       number = 6,
                                      repeats = 6, classProbs = TRUE)) # todo check these args

fit <- classify(data=data.trainS4, method="svmRadial",
                preProcessing="deseq-vst", ref="male",
                control = trainControl(method = "repeatedcv", 
                                       number = 6,
                                       repeats = 6)) # todo check these args


show(fit)
pred.svm <- predict(fit, data.testS4)
pred.svm <- relevel(pred.svm, ref = "male")
actual <- relevel(classts$condition, ref = "male")
tbl <- table(Predicted = pred.svm, Actual = actual)
confusionMatrix(tbl, positive = "male")
trained(fit)
plot(fit)


# ----- just use their preprocessing setup ---- #
rawData <- data
transformedData <- NULL

preProcessing <- match.arg(preProcessing)
preProcessData <- function(data, preProcessing){
  if (preProcessing == "deseq-vst"){
    normalization <- "deseq"
    transformation <- "vst"
    
    dataSF <- estimateSizeFactors(data)    # Estimate Size factors:
    dataDisp <- estimateDispersions(dataSF, fitType = "local")  # Estimate dispersions:
    transformedData <- varianceStabilizingTransformation(dataDisp, fitType = "local")
    
    dataVST <- t(as.matrix(assay(transformedData)))
    input <- dataVST   ## Input data from transformed expression data.
    output <- colData(data)[ ,class.labels]   ## Output: class labels of samples.
    
    trainParameters <- list(disperFunc = dispersionFunction(dataDisp),
                            sizeFactors = sizeFactors(dataDisp))
    
  } else if (preProcessing == "deseq-rlog"){
    normalization <- "deseq"
    transformation <- "rlog"
    
    dataSF <- estimateSizeFactors(data)    # Estimate Size factors:
    dataDisp <- estimateDispersions(dataSF, fitType = "local")  # Estimate dispersions:
    
    transformedData = rlog(dataDisp, blind = FALSE)  ### RLOG donusumu uygulanmis olan bu nesnenin de MLSeq objesi icerisinde predictClassify'a gonderilmesi gerekiyor.
    dataRLOG = t(assay(transformedData))
    input <- dataRLOG   ## Input data from transformed expression data.
    output <- colData(data)[ ,class.labels]   ## Output: class labels of samples.
    
    # dataexpTrainRLOG = t(assay(trS4expRLOG))
    # betaPriorVar = attr(trS4expRLOG,"betaPriorVar")
    # intercept = mcols(trS4expRLOG)$rlogIntercept
    
    trainParameters <- list(betaPriorVar = attr(dataRLOG,"betaPriorVar"),
                            intercept = attr(dataRLOG, "intercept"),
                            dispFit = mcols(transformedData)$dispFit)
    
  } else if (preProcessing == "deseq-logcpm"){
    normalization <- "deseq"
    transformation <- "logcpm"
    
    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "RLE")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.
    
    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes
    
    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes
    
    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]   ## Output: class labels of samples.
    
    trainParameters <- list(NULL)
    
  } else if (preProcessing == "tmm-logcpm"){
    normalization <- "tmm"
    transformation <- "logcpm"
    
    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.
    
    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes
    
    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes
    
    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]   ## Output: class labels of samples.
    
    # This function is used to find reference sample.
    # Codes are copied from edgeR and modified here.
    findRefSample <- function (rawCounts, lib.size, p = 0.75){
      y <- t(t(rawCounts) / lib.size)
      f75 <- apply(y, 2, function(x){
        quantile(x, p = p)
      })
      refSample <- which.min(abs(f75-mean(f75)))
      return(refSample)
    }
    
    refSampleID <- findRefSample(rawCounts, lib.size = colSums(rawCounts), p = 0.75)
    trainParameters <- list(refSample = rawCounts[ ,refSampleID])
    
  } else if (preProcessing == "logcpm"){
    normalization <- "none"
    transformation <- "logcpm"
    
    rawCounts = counts(data, normalized = FALSE)
    countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
    countsDGE.normalized <- calcNormFactors(countsDGE, method = "none")   ## RLE: DESeq mantigi ile normalize ediyor.
    countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.
    
    rawData <- countsDGE
    classes <- colData(data)[ ,class.labels]
    rawData[[class.labels]] <- classes
    
    transformedData <- list(normalizedData = countsDGE.normalized, transformedData = countsDGE.transformed)
    transformedData[[class.labels]] <- classes
    
    input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
    output <- rawData[[class.labels]]  ## Output: class labels of samples.
    
    trainParameters <- list(NULL)
  }
}


setupPredict <- function(model, test.data){
  test.pred <- NULL
  
  ## Prediction steps for "caret" based classifications.
    if (normalization(model) == "deseq" & transformation(model) == "vst"){  # "deseq-vst"
      counts.test <- counts(test.data)
      counts.train <- counts(model@inputObject$rawData)    ### Raw Counts icin getter yazilinca burada kullanilacak.  rawCounts(...) adinda bir getter kullanilabilir.
      
      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = counts.test / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))
      
      test.dataSF <- estimateSizeFactors(test.data)    # Estimate Size factors:
      sizeFactors(test.dataSF) <- sizeF.ts             # Replace size factors with size factors which are estimates using train set parameters.
      
      ## Change dispersion function of test data with dispersion function of train data
      dispersionFunction(test.dataSF) <- trainParameters(model)$disperFunc
      
      transformedData <- varianceStabilizingTransformation(test.dataSF, fitType = "local", blind = FALSE)
      
      dataVST <- t(as.matrix(assay(transformedData)))
      input <- dataVST   ## Input data from transformed expression data.
      
    } else if (normalization(model) == "deseq" & transformation(model) == "rlog"){   # "deseq-rlog"
      counts.test <- counts(test.data)
      counts.train <- counts(model@inputObject$rawData)    ### Raw Counts icin setter yazilinca burada kullanilacak.  rawCounts(...)
      
      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = counts.test / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))
      
      test.dataSF <- estimateSizeFactors(test.data)    # Estimate Size factors:
      sizeFactors(test.dataSF) <- sizeF.ts
      
      test.dataDisp <- estimateDispersions(test.dataSF, fitType = "local")
      
      ## Required train parameters for test data
      mcols(test.dataDisp)$dispFit <- trainParameters(model)$dispFit
      betaPrior <- trainParameters(model)$betaPrior
      intercept <- trainParameters(model)$intercept
      
      test.dataRLOG <- rlog(test.dataDisp, blind = FALSE,
                            intercept = intercept, betaPriorVar = betaPrior)
      
      dataRLOG = t(assay(test.dataRLOG))
      input <- dataRLOG   ## Input data from transformed expression data.
      
    } else if (normalization(model) == "deseq" & transformation(model) == "logcpm"){  ## deseq-logcpm
      
      rawCounts = counts(test.data, normalized = FALSE)
      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
      counts.train <- counts(model@metaData$metaData@rawData.DESeqDataSet)    ### Raw Counts icin setter yazilinca burada kullanilacak.  rawCounts(...)
      
      #Calculation of test set size factors using geometric means from training data
      #Genes in the row, samples in the column
      geomts = rawCounts / exp(rowMeans(log(counts.train)))   ## Geometric mean of test data using estimators from train data
      sizeF.ts = apply(geomts, 2, function(x) median(x))
      
      RLE <- sizeF.ts / colSums(rawCounts)
      normFactors <- RLE / exp(mean(log(RLE)))
      
      countsDGE.normalized <- calcNormFactors(countsDGE, method = "RLE")   ## RLE: DESeq mantigi ile normalize ediyor.
      countsDGE.normalized$samples$norm.factors <- normFactors   ### Replace normFactors using normFactors obtained from train set parameters.
      countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek
      
      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
      
    } else if (normalization(model) == "none" & transformation(model) == "logcpm"){
      rawCounts = counts(test.data, normalized = FALSE)
      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
      countsDGE.transformed <- cpm(countsDGE, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek
      
      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
      
    } else if (normalization(model) == "tmm" & transformation(model) == "logcpm"){
      
      rawCounts = counts(test.data, normalized = FALSE)
      referenceSample <- trainParameters(model)$refSample
      
      ## Chech if the feature names of reference sample are in the same order with rawCounts.
      if (identical(rownames(rawCounts), names(referenceSample))){
        rawCounts <- cbind(rawCounts, referenceSample)
      } else {
        stop(warning("Reference sample either does not have same features or the features are not in the same order as test set. Calculation stops.."))
      }
      
      countsDGE <- DGEList(counts = rawCounts, genes = rownames(rawCounts))
      
      ## Reference column is selected from train set.
      countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM", refColumn = ncol(countsDGE))   ## RLE: DESeq mantigi ile normalize ediyor.
      countsDGE.transformed <- cpm(countsDGE.normalized, log = TRUE, prior.count = 1)   ### prior Count daha sonra duzenlenecek.
      
      input <- t(countsDGE.transformed)   ## Input data from transformed expression data.
      input <- input[-nrow(input), ]   ## Remove reference sample from test data.
    }
    
    test.pred = predict(trained(model), input, ...)    ## Predicted classes.
    
  

  
  # res <- list(MLSeqObject = model, Predictions = test.pred)
  return(test.pred)

}



# TODO - can I set this up to do a wider range of values?

# voomNSC
#  - voom transformation+ NSC method
#  for these normalize with deseq, and then run model

# types of classifiers
# - continuous (caret) - trainControl
# - discrete - discreteControl
# - voom-based - voomControl

# ----
# "trimmed mean of M values" (TMM)
#   - trims the data in upper log-fold changes to minimize btw-sample / intensity changes

# log-cpm
#


# sparse classifiers for determining possible biomarkers
# "selectedGenes()" <-- how you get the list