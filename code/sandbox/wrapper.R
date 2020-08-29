
  

    clusters <- lapply(AllGSEmetadata(), function(X){
                     return(GEOclustering(X))
                   })

    Predictions <- lapply(names(GSEclusters()), function(X){LabelClusters(GSEclusters()[[X]], AllGSEmetadata()[[X]], model())})
      names(Predictions) <- names(GSEclusters())

       FixedPredictions <- mapply(FixLabels, GSEclusters(), clusterLabels(), input$analysis_type, SIMPLIFY = F)
      
      invalid <- unlist(lapply(FixedPredictions, function(X){return(X$confidence[1] == "Invalid")}))
    FixedLabels <- FixedPredictions[!invalid])
    MatchedPairs <- lapply(names(FixedLabels), function(X){
                     res <- MatchClusters(FixedLabels[[X]], AllGSEmetadata()[[X]])
                     names(res) <- unlist(lapply(res, function(X){paste(X$Mut, X$WT, collapse="_")}))
                     return(res)
                   })
                   names(MatchedPairs) <- names(FixedLabels)
  
  

 MPCs <- MatchedPairs[unlist(lapply(MatchedPairs, length)) > 0])

###########
# next is FPMCs

  MatchedPairs <- eventReactive(GSEclusters(), {
    
    print("Matching clusters")
    
    
    clusterLabels <- reactive({
      Predictions <- lapply(names(GSEclusters()), function(X){LabelClusters(GSEclusters()[[X]], AllGSEmetadata()[[X]], model())})
      names(Predictions) <- names(GSEclusters())
      return(Predictions)
    })
    
    FixedLabels <- reactive({
      
      FixedPredictions <- mapply(FixLabels, GSEclusters(), clusterLabels(), input$analysis_type, SIMPLIFY = F)
      
      invalid <- unlist(lapply(FixedPredictions, function(X){return(X$confidence[1] == "Invalid")}))
      return(FixedPredictions[!invalid])
    })
    
    withProgress(message = "IN PROGRESS:",
                 detail = "Matching Clusters: Initialising", value = 0.1, {
                   Pairs <- lapply(names(FixedLabels()), function(X){
                     incProgress(amount = 0.8/length(FixedLabels()), detail = paste0("Matching Clusters: ", X))
                     res <- MatchClusters(FixedLabels()[[X]], AllGSEmetadata()[[X]])
                     names(res) <- unlist(lapply(res, function(X){paste(X$Mut, X$WT, collapse="_")}))
                     return(res)
                   })
                   names(Pairs) <- names(FixedLabels())
                 })
    return(Pairs)
  })
  
  GSEplatforms <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      gpls <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$gpl
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  GSEspecies <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      species <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$organism_ch1
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  
  GSEchannels <- reactive({
    res <- lapply(names(MatchedPairs()), function(x){
      X <- AllGSEmetadata()[[x]]
      
      channels <- unique(unlist(lapply(X$GSMMeta, function(x){
        x$label_ch2
      })))
    })
    names(res) <- names(MatchedPairs())
    return(res)
  })
  
  
  output$numGSEs <- renderText(
    paste0('You have ', length(MatchedPairs()), ' GSE that match your selected experiment type.')
  )
  
  