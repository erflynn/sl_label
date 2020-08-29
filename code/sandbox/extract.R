# metadata functions

############

# This function grabs all the metadata related to a GSE (experiment) ID, includng the associated GSM (sample) IDs
grabAllMetadataGSE <- function(GSEid) {
  # Fetch the GSE metadata and make it a list
  GSEMeta <- as.list(dbGetQuery(metadata,paste("select * from gse where gse.gse =\"", GSEid,"\"", "\n", sep="")))
  # Get the GSM IDs contained within this GSE
  GSEtoGSM<- grabGSMids(GSEid)
  # Get the GSM metadata
  GSMMeta <- lapply(GSEtoGSM, grabAllMetadataGSM)
  # Add the GSM metadata to the GSE metadata
  if(length(unlist(GSMMeta)) == 0){return(NA)}
  
  GSEMeta$GSMMeta <- GSMMeta
  # Return the GSE metadata list
  return(GSEMeta)
}

# This function grabs all metadata related to a GSM (sample) ID
grabAllMetadataGSM <- function (gsmid){
  dbGetQuery(metadata, paste("select * from gsm where gsm.gsm =\"", gsmid,"\"\n", sep=""))
}

# This function finds the GSM (sample) IDs associated with a GSE (experiment) ID 
grabGSMids <- function(GSEid){
  gsm <- unlist(dbGetQuery(metadata,paste("select gse_gsm.gsm from gse_gsm where gse_gsm.gse =\"", GSEid,"\"", "\n", sep="")))
  return(gsm)
}

# This function finds the GPL (platform) ID associated with a GSM (sample) ID
grabGPL <- function (gsmid){
  gplid <- unlist(dbGetQuery(metadata,paste("select gsm.gpl from gsm where gsm.gsm =\"", gsmid,"\"", "\n", sep="")))
  return(gplid)
}

## select metadata from multiple GSE ids, returns a data.frame instead of a list
grabMultipleGSE <- function (gseids){ 
  dbGetQuery(metadata, paste("select * from gse where gse IN (\'" ,paste(gseids, collapse="\',\'"), "\')", sep="")) 
}

## select metadata from multiple GSM ids, returns a data.frame instead of a list
grabMultipleGSM <- function (gsmids){
  dbGetQuery(metadata, paste("select * from gsm where gsm.gsm IN (\'" ,paste(gsmids, collapse="\',\'"), "\')", sep="")) 
}


################################################################################################################################################
# GEORACLE FUNCTIONS
################################################################################################################################################



get.title.table <- function(GSM){
  title <- GSM$title
  o.title <- title
  if(!grepl("(wt)|(ko)(_)?\\d*?$",title, ignore.case = T) & !grepl("((\\s|\\_|\\-)[0-9]{1,2}(e|h|w|d|p|(pd)|(day)|(week)|(ed)|(hr)|(hour))$)|((\\s\\_\\-)(e|h|w|d|p|(pd)|(day)|(week)|(ed)|(hr)|(hour))[0-9]{1,2}$)", title, ignore.case = T)){ 
    title <- unlist(strsplit(gsub("[!a-z!0-9!A-Z][0-9a-zA-Z]{,2}$","",title), "\\,|\\_"))
  }
  if(identical(title, character(0))){ title <- o.title}
  title <- unlist(strsplit(gsub("\\s|\\-","_",title), "\\,|\\_"))
  title <- title[!title == ""]
  return(title)
}

get.characteristics.table <- function(GSM){
  characteristics <- unlist(strsplit(GSM$characteristics_ch1, ";\t|\\, "))
  characteristics.split <- do.call(rbind, lapply(characteristics, function(X){return(unlist(strsplit(X, ":\\s?")))}))
  
  if(as.numeric(ncol(characteristics.split)) > 1) {
    characteristics.table <- as.data.frame(t(characteristics.split[,2]))
    colnames(characteristics.table) <- characteristics.split[,1]
  } else {
    characteristics.table <- as.data.frame(characteristics.split)
    if(nrow(characteristics.table)>1) {
      characteristics.table <- as.data.frame(t(characteristics.table))
    }
  }
  return(characteristics.table)
}

#Calculates combined table
get.combined.table <- function(GSM) {
  characteristicsTable <- get.characteristics.table(GSM)
  titleTable <- get.title.table(GSM)
  combined.table <- cbind(t(titleTable),characteristicsTable)
  return(combined.table)
}

###########################################################
###########################################################


get.source <- function(GSM) {
  source <- GSM$source_name_ch1
  return(source)
}


## takes in GSM metadata
get.combined.table2 <- function(GSM) {
  characteristicsTable <- get.characteristics.table(GSM)
  titleTable <- get.title.table(GSM)
  souceInfo <- get.source(GSM)
  combined.table <- cbind(t(titleTable),characteristicsTable,souceInfo)
  return(combined.table)
}

get.min.diss <- function(Table){
  
  Table <- data.matrix(Table)
  Table[is.na(Table)] <- "a"
  Table <- data.frame(Table)
  diss.mat <- data.matrix(daisy(Table))
  
  min.dis <- lapply(1:nrow(diss.mat), function(X){
    return(unlist(which(diss.mat[X,] == min(diss.mat[X,]))))
  })
  return(min.dis)
}

TitleClustering <- function(gseMetadata) {
  
  print(paste0("Clustering ", gseMetadata$gse))
  TitleTable <- data.frame(do.call(rbind, lapply(gseMetadata$GSMMeta, get.title.table)))
  
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  if(sum(duplicated(GSMtitles)) == 0){
    rownames(TitleTable) <- GSMtitles
  } else {
    rownames(TitleTable) <- names(GSMtitles)
  }
  # finding least dissimilarity = best match
  title.min.dis <- get.min.diss(TitleTable)
  
  tNamesInClusters <- lapply(unique(title.min.dis), function(X){return(GSMtitles[X])})
  
  return(tNamesInClusters)  
}

#Group GSMs based on characteristics
CharacteristicsClustering <- function(gseMetadata) {
  
  CharTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.characteristics.table)))
  
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  if(sum(duplicated(GSMtitles)) == 0){
    rownames(CharTable) <- GSMtitles
  } else {
    rownames(CharTable) <- names(GSMtitles)
  }# finding least dissimilarity = best match
  char.min.dis <- get.min.diss(CharTable)
  
  cNamesInClusters <- lapply(unique(char.min.dis), function(X){return(GSMtitles[X])})
  return(cNamesInClusters)  
}

#Group GSMs based on both title and characteristics
CombinedClustering <- function(gseMetadata) {
  combinedTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.combined.table2)))
  GSMtitles <- unlist(lapply(gseMetadata$GSMMeta, function(X){return(X$title)}))
  names(GSMtitles) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  combined.min.dis <- get.min.diss(combinedTable)
  fNamesInClusters <- lapply(unique(combined.min.dis), function(X){return(GSMtitles[X])})
  return(fNamesInClusters)  
}

#Check if clustering is valid 
informative.clustering <- function(CL){
  if(length(CL) > 1 & length(CL) <  (length(unlist(CL)) - 1) ){
    return("TRUE")
  }
  return("FALSE")
}

#Check for invalid titles which are purely numeric
non.numeric.title <- function(GSM) {
  title <- GSM$title
  title <- gsub("[!a-z!0-9!A-Z][0-9a-zA-Z]{,2}$","",title)
  title <- gsub("\\s|\\-","",title)
  if(grepl("^[0-9]+$", title)) {
    return("FALSE")
  }
  return("TRUE")
}

#Decide whether to use title clustering or characteristics clustering 
GEOclustering <- function(gseMetadata){
  
  title.res <- TitleClustering(gseMetadata)
  char.res <- CharacteristicsClustering(gseMetadata)
  combined.res <- CombinedClustering(gseMetadata)
  
  names(title.res) <- unlist(lapply(title.res, function(X){ return( paste(names(X), collapse="_"))}))
  names(char.res) <- unlist(lapply(char.res, function(X){ return( paste(names(X), collapse="_"))}))
  names(combined.res) <- unlist(lapply(combined.res, function(X){ return( paste(names(X), collapse="_"))}))
  
  tinfo <- informative.clustering(title.res)
  non.numeric <- unlist(lapply(gseMetadata$GSMMeta, non.numeric.title))
  if("FALSE" %in% non.numeric) {
    tvalid <- "FALSE"
  } else {
    tvalid <- "TRUE"
  }
  
  cinfo <- informative.clustering(char.res)
  
  if (tinfo == "TRUE" && cinfo == "TRUE") {
    #Check whether title and characteristics clusterings match
    factor.char.res <- lapply(char.res, as.factor)
    level.char.res <- lapply(factor.char.res, function(X){
      combineLevels(X,c(X[1:length(X)]), newLabel = toString(X[1]))})
    #Unlist cluster 
    unlisted.char <- as.factor(unlist(level.char.res))
    tClustersInCharClusters <- lapply(title.res, function(X){
      return(unique(unlisted.char[match(X, as.factor(unlist(char.res)))]))
    })
    
    maxNumberOfCharPerT <- max(unlist(lapply(tClustersInCharClusters, length)))
    
    maxNumberOfTPerChar <- max(table(unlist(tClustersInCharClusters)))
    
    
    if(maxNumberOfCharPerT == 1 && maxNumberOfTPerChar == 1) {
      #Matching title and characteristics clustering
      return(title.res)
    } else if(tvalid == "FALSE") {
      #Numeric titles 
      return(char.res)
    } else {
      #Title and characteristics clustering not matching
      factor.combined.res <- lapply(combined.res, as.factor)
      level.combined.res <- lapply(factor.combined.res, function(X){
        combineLevels(X,c(X[1:length(X)]), newLabel = toString(X[1]))})
      
      unlisted.combined <- as.factor(unlist(level.combined.res))
      tClustersInCombinedClusters <- lapply(title.res, function(X){
        return(unique(unlisted.combined[match(X, as.factor(unlist(combined.res)))]))
      })
      
      charClustersInCombinedClusters <- lapply(char.res, function(X){
        return(unique(unlisted.combined[match(X, as.factor(unlist(combined.res)))]))
      })
      
      maxNumberOfCombinedPerT <- max(unlist(lapply(tClustersInCombinedClusters, length)))
      
      maxNumberOfTPerCombined <- max(table(unlist(tClustersInCombinedClusters)))
      
      maxNumberOfCombinedPerChar <- max(unlist(lapply(charClustersInCombinedClusters, length)))
      
      maxNumberOfCharPerCombined <- max(table(unlist(charClustersInCombinedClusters)))
      
      if((maxNumberOfCombinedPerT == 1 && maxNumberOfTPerCombined == 1)|
         (maxNumberOfCombinedPerChar == 1 && maxNumberOfCharPerCombined == 1)){
        return(combined.res)
      } else {
        return(combined.res)
      }
    }
  } else if(tinfo == "TRUE" && cinfo == "FALSE") {
    #Invalid characteristics clustering
    return(title.res)
  } else if(tinfo == "FALSE" && cinfo == "TRUE") {
    #Invalid title clustering 
    return(char.res)
  } else {
    return("Discard - Neither informative")
  }
}

###########################################################
###########################################################

#Extract title features for GSM
get.title.features <- function(cluster) {
  title <- paste(cluster, collapse = " ")
  clusterTitleFeatures <- levels(factor(unlist(strsplit(gsub("\\,|:|\\.|;|\\_|\\)|\\("," ",title), "\\s{1,3}"))))
  return(clusterTitleFeatures)
}

#Gets characteristics features for GSE
get.characteristics.features <- function(GSEclusters, CurGSEmetadata) {
  #Get clusters in GSM codes
  clusterCode <- lapply(GSEclusters,names)
  
  GSMChar <- lapply(lapply(clusterCode, function(x){
    lapply(x, function(y){
      lapply(CurGSEmetadata$GSMMeta, function(z) {
        z$characteristics_ch1[z$gsm == y]
      })
    })
  }), unlist)
  
  clusterChar <- lapply(GSMChar, function(x){ 
    paste(x, collapse = ";")
  })
  
  #Tokenise keywords of characteristics 
  clusterCharFeatures <- lapply(clusterChar, function(x) {
    x2 <- strsplit(x, ";")[[1]]
    x3 <- gsub(".*:(.*)","\\1",x2)
    levels(factor(unlist(strsplit(gsub("\\,|:|\\.|;|\\_|\\)|\\("," ",x3), "\\s+"))))
  })
  
  return(clusterCharFeatures)
}

#Return list with info on cluster GSMS, title keywords and characteristic keywords 
combine.cluster.features <- function(GSEclusters, CurGSEmetadata) {
  
  #Get features and pu into string
  clusterChar <- lapply(get.characteristics.features(GSEclusters, CurGSEmetadata), function(x){ 
    paste(x, collapse = " ")
  })
  
  clusterTitle <- lapply(lapply(GSEclusters, get.title.features), function(x){ 
    paste(x, collapse = " ")
  })
  
  clusterGSM <- lapply(lapply(GSEclusters, names), function(x){ 
    paste(x, collapse = " ")
  })
  
  combinedTable <- cbind(clusterGSM, clusterTitle, clusterChar)
  
  return(combinedTable)
}

############################################

#Load control and treatment features
get.features <- function(txtFile) {
  featureTxt <- read.table(txtFile)
  features <- as.vector(unlist(featureTxt))
}

############################################

#Check presence of features
checkPresence <- function(GSEFeatures) {
  
  combined.features <- get.features(Features.file)
  
  feat.matrix <- do.call(rbind, lapply(GSEFeatures, function(cluster){
    feat.vector <- unlist(lapply(combined.features, function(feat){
      grepl(feat, cluster, ignore.case = T)
    }))
  }))
  
  colnames(feat.matrix) <- combined.features
  return(feat.matrix)
}

#Produce logic matrix to inidicate presence of features in title and characteristics
get.logic <- function(GSEclusters, CurGSEmetadata) {
  
  GSEFeatures <- combine.cluster.features(GSEclusters, CurGSEmetadata)
  
  #Check presence of features in title
  titleFeatures <- checkPresence(GSEFeatures[,2])
  
  #Check presence of features in characteristics
  charFeatures <- checkPresence(GSEFeatures[,3])
  
  titleCharMatrix <- cbind(titleFeatures, charFeatures)
  return(titleCharMatrix)
}

###############


#Obtain matrix
get.svm.m <- function(GSEclusters, CurGSEmetadata) {
  
  m <- get.logic(GSEclusters, CurGSEmetadata)
  
  #Edit column names
  colnames(m) <- gsub("/","neg",colnames(m))
  colnames(m) <- gsub("-","\\.",colnames(m))
  colnames(m) <- gsub("\\?|\\|","\\.",colnames(m))
  
  #Distinguish between characteristics and title features
  colnames(m)[(ncol(m)/2+1):ncol(m)] <- paste(colnames(m)[(ncol(m)/2+1):ncol(m)],".char", sep='')
  
  return(m)
}

#Obtain predictions through SVM
get.preds <- function(data, model) {
  
  #Define dependent variable
  d <- data.frame(data)
  
  predictions <- predict(model, d, probability = T)
  return(predictions)
}




LabelClusters <- function(GSEclusters, CurGSEmetadata, model){
  data <- get.svm.m(GSEclusters, CurGSEmetadata)
  prediction <- get.preds(data, model)
  
  return(prediction)
}

#Reverse label
reverse.label <- function(oldLabel) {
  newLabel <- c("FALSE", "TRUE")[c("FALSE", "TRUE")!= as.vector(oldLabel)]
  return(newLabel)
}

FixLabels <- function(GSEclusters, predictions, analysis_type = "All experiments"){
  
  
  GSEClassInfo <- data.frame(predLabel = predictions, probs = attr(predictions, "probabilities"), row.names = lapply(lapply(GSEclusters, names), paste, collapse=" "))
  
  reverseRow <- vector()
  assessClass <- vector()
  threshold <- 0.83
  
  #Check if number of classes
  if(length(unique(GSEClassInfo$predLabel)) == 1) {
    #If one class:
    
    #Select probability values to work with 
    if(unique(GSEClassInfo$predLabel) == "FALSE") {
      GSEClassInfo$probability <- GSEClassInfo$probs.FALSE
    } else{
      GSEClassInfo$probability <- GSEClassInfo$probs.TRUE
    }
    
    #Check if all probability values equal 
    if (!all(GSEClassInfo$probability == GSEClassInfo$probability[1])) {
      
      #Check for probability values above threshold
      highConf <- rownames(GSEClassInfo)[GSEClassInfo$probability >= threshold]
      
      #Label high confidence clusters
      if(length(highConf) > 0 & length(highConf) < dim(GSEClassInfo)[1]) {
        for(i in 1:length(highConf)){
          assessClass[i] <- "High"
          names(assessClass)[i] <- highConf[i]
        }
      } else if (length(highConf) == nrow(GSEClassInfo)) {
        #Clear high confidence vector if all cluster probabilities above threshold
        length(highConf) <- 0
      }
      
      #Reverse labels with lowest confidence 
      reverseRow <- rownames(GSEClassInfo)[which(GSEClassInfo$probability == min(GSEClassInfo$probability))]
      GSEClassInfo[reverseRow,]$predLabel <- sapply(GSEClassInfo[reverseRow,]$predLabel,reverse.label)
      
      for(i in 1:length(reverseRow)){
        assessClass[i+length(highConf)] <- "Low"
        names(assessClass)[i+length(highConf)] <- reverseRow[i]
      }
      
      #Label everything neither high nor low as medium
      mediumClusters <- rownames(GSEClassInfo)[which(!rownames(GSEClassInfo) %in% names(assessClass))]
      j <- length(assessClass)
      if (length(mediumClusters) != 0) {
        for(i in 1:length(mediumClusters)){
          assessClass[i+j] <- "Medium"
          names(assessClass)[i+j] <- mediumClusters[i]
        }
      }
    } else {
# actually, lets always choose one to randomly be the "control" if we cant detect it...     
#      if(analysis_type == "All experiments"){
      
          ## if we want to include all GSE, we will randomly permute one cluster to be different
        assessClass <- rep("Low", nrow(GSEClassInfo))
        randomCluster <- sample(1:nrow(GSEClassInfo), 1)
        GSEClassInfo$predLabel[randomCluster] <- setdiff(c(T,F),GSEClassInfo$predLabel[randomCluster])
        
#      }else{
        #All probabilities equal - invalid clusters 
#        assessClass <- rep("Invalid", nrow(GSEClassInfo))
#      }
      names(assessClass) <- rownames(GSEClassInfo)
    }
  } else {
    #If both classes present:
    #Check for clusters that are high confidence
    highConfMut <- rownames(GSEClassInfo)[GSEClassInfo$probs.FALSE >= threshold]
    highConfWT <- rownames(GSEClassInfo)[GSEClassInfo$probs.TRUE >= threshold]
    highConf <- c(highConfMut,highConfWT)
    
    #Label high confidence clusters
    if(length(highConf) > 0) {
      for(i in 1:length(highConf)){
        assessClass[i] <- "High"
        names(assessClass)[i] <- highConf[i]
      }
    } 
    
    #Label medium confidence clusters  
    mediumClusters <- rownames(GSEClassInfo)[which(!rownames(GSEClassInfo) %in% names(assessClass))]
    j <- length(assessClass)
    if (length(mediumClusters) != 0) {
      for(i in 1:length(mediumClusters)){
        assessClass[i+j] <- "Medium"
        names(assessClass)[i+j] <- mediumClusters[i]
      }
    }
  }
  
  GSEClassInfo$confidence = assessClass[match(rownames(GSEClassInfo), names(assessClass))]
  
  return(GSEClassInfo)
}



###################################################################
#Get GSM clusters and labels
get.cluster.labels <- function(GSEId,class.res) {
  info.table <- class.res[grepl(GSEId,rownames(class.res)),][,1:2]
  colnames(info.table) <- c("GSM", "Label")
  
  return(info.table)
}

#check if multi-control analysis required
check.multi <- function(clusterGSMLabel) {
  if((nrow(clusterGSMLabel) == 2) | (sum(clusterGSMLabel[,2] == "FALSE") == 1)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
######################################################
#Get label vector describing class or group
get.labels <- function(labelledClusters) {
  splitGSM <- lapply(labelledClusters,function(x){unlist(strsplit(as.character(x)," "))})
  GSMLabels <- sub("^([[:alpha:]]*).*", "\\1", names(unlist(splitGSM)))
  names(GSMLabels) <- unlist(splitGSM) 
  return(GSMLabels)
}

#Get table with title, characteristics and source keywords plus class and subgroup labels
get.label.table <- function(clusterGSMLabel,gseMetadata) {
  #Get combined table
  combinedTable <- data.frame(do.call(rbind.fill, lapply(gseMetadata$GSMMeta, get.combined.table2)))
  rownames(combinedTable) <- do.call(rbind, lapply(gseMetadata$GSMMeta, function(x) {x$gsm}))
  
  #Generating class vector for GSMs
  mC <- as.data.frame(clusterGSMLabel)
  classSplit <- split(mC$GSM,mC$Label)
  classLabel <- get.labels(classSplit)
  
  #Generating subgroup vector for GSMs
  ms <- clusterGSMLabel[,1]
  names(ms) <- letters[1:length(ms)]
  subgroupLabel <- get.labels(ms)
  
  labelT <- cbind(combinedTable,ClassLabels = classLabel[rownames(combinedTable)], SubgroupLabels = subgroupLabel[rownames(combinedTable)])
  
  return(labelT)
}

#Get dissimilarity matrix of all GSMs
get.diss.mat <- function(Table){ #Put in combinedTable
  Table <- data.matrix(Table)
  Table[is.na(Table)] <- "a"
  Table <- data.frame(Table)
  diss.m <- data.matrix(daisy(Table))
  return(diss.m)
}

#Input cluster of GSM titles as string and output list of matched pair
pair.clusters <- function(GSMCluster,labelTable, CurGSEmetadata) {
  clusterGSM <- strsplit(GSMCluster," ")[[1]]
  clusterLabel <- unique(labelTable[GSMCluster,]$predLabel)
  if(clusterLabel == "FALSE") {
    return(NA)
  }
  
  combinedTable <- data.frame(do.call(rbind.fill, lapply(CurGSEmetadata$GSMMeta, get.combined.table2)))
  rownames(combinedTable) <- do.call(rbind, lapply(CurGSEmetadata$GSMMeta, function(x) {x$gsm}))
  #Obtain dissimilarity matrix of cluster with potential GSM matches 
  diss.mat <- get.diss.mat(combinedTable)
  
  diss.mat.cluster <- diss.mat[rownames(diss.mat) %in% clusterGSM,, drop=FALSE]
  
  
  #Obtain GSM of potential matches
  pMatch <- rownames(labelTable)[labelTable$predLabel != clusterLabel]
  names(pMatch) <- pMatch
  pDisss <- lapply(pMatch, function(CL){
    GSMs <- unlist(strsplit(CL, " "))
    return(diss.mat.cluster[,GSMs, drop = FALSE])
  })
  
  pMatchMeanDiss <- lapply(pDisss, mean)
  
  matchedClusters <- pMatchMeanDiss[which(pMatchMeanDiss == min(as.numeric(pMatchMeanDiss)))]
  
  ##
  ## Currently in the case of a tie just return the first one
  matchedGSMs <- names(matchedClusters[1])
  
  #Assess confidence 
  if(length(matchedClusters) == 1) {
    confidence <- "High"
  } else {
    confidence <- "Low"
  }    
  
  result <- list(GSMCluster,matchedGSMs,confidence)
  
  names(result) <- c("Mut","WT","Confidence")
  
  #Return paired labelled clusters 
  return(result)
}

#Match clusters for simple GSEs
simple.pair <- function(GSMCluster,labelTable) {
  clusterLabel <- unique(labelTable[GSMCluster,]$predLabel)
  if(clusterLabel == "FALSE") {
    return(NA)
  }
  
  matchedGSMs <- rownames(labelTable)[labelTable$predLabel == "FALSE"]
  
  result <- list(GSMCluster,matchedGSMs,"High")
  names(result) <- c("Mut","WT","Confidence")
  #Return paired labelled clusters 
  return(result)
}

#save matched GSMs in text file 
save.match <- function(GSMres,GSEId) {
  fName <- paste(GSMres[["Mut"]],collapse = "-")
  res1 <- lapply(GSMres, function(x){
    paste(x,collapse = " ")
  })
  writeData <- as.data.frame(do.call(rbind,res1))
  write.table(writeData, file = paste(GSEId,"-",fName,".txt",sep =""),col.names = FALSE)
}

#Get match and confidence
get.match.outcome2 <- function(labelTable, CurGSEmetadata) { #Pass through matrix with labelled clusters 
  multi <- check.multi(labelTable)
  
  if(multi) {
    outcome <- lapply(rownames(labelTable[labelTable$predLabel == "TRUE",,drop = FALSE]),function(x) {
      pair.clusters(x,labelTable, CurGSEmetadata)
    })
  } else {
    outcome <- lapply(rownames(labelTable[labelTable$predLabel == "TRUE",,drop = FALSE]),function(x) {
      simple.pair(x,labelTable)
    })
  } 
  return(outcome)
}




MatchClusters <- function(clusterLabels, CurGSEmetadata){
  get.match.outcome2(clusterLabels, CurGSEmetadata)
}

