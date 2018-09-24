# Meta (GET)
# This is a Shiny web application for meta-analysis in SyRF
# by Qianying Wang @CAMARADES
# server.R
# https://qianying.shinyapps.io/Meta/

library(shiny)
library(metafor)
library(meta)
library(shinythemes)
library(dplyr)
library(plotly)
library(RCurl)
library(ggplot2)
library(colourpicker)


shinyServer(function(input, output, session) {
  
  # dataset <- reactive({
  #   query <- parseQueryString(session$clientData$url_search)
  # 
  #   if (is.null(query[["csvpath"]])) { stop() }
  #   link <- query[["csvpath"]]
  #   # link <- "http://www.dcn.ed.ac.uk/camarades/Behaviour outcomeData.csv"
  #   # myfile <- getURL(link, ssl.verifyhost=F, ssl.verifypeer=F)
  #   dataset <- read.csv(link, header=T)
  #   
  #    # dataset <- read.csv("http://www.dcn.ed.ac.uk/camarades/Behaviour outcomeData.csv", header = T)
  # })
  dataset <- reactive({
    inFile = input$file
    if (is.null(inFile)) { stop() }
    dataset <- read.csv(inFile$datapath, header = T, row.names = NULL)
    # dataset <- na.omit(dataset)
    # dataset <- read.csv("D:/shinyapp/datafile/Infarct.csv", header=T)
  })
  
  # ------------------ Record: yi & vi --------------------
  record <- reactive({
    
    dat <- dataset()
    dat[dat[,"GreaterIsWorse"]=="TRUE", "GreaterIsWorse"]=1
    dat[dat[,"GreaterIsWorse"]==0, "GreaterIsWorse"]=-1
    
    fix = 18 # number of fixed columns in dataset                   
    
    record = c()
    
    type = input$CompareType
    mea = input$EffectMeasure
    
    
    ExpList = split(dat, dat[,"ExperimentID"], drop = T)
    
    for (i in 1:length(ExpList)){ # within i-th experiment
      
      samp = ExpList[[i]] 
      OutList = split(samp, samp[,"OutcomeId"], drop = T)
      
      for (j in 1:length(OutList)){ # within j-th outcome
        
        out = OutList[[j]] 
        stem = out[(out[,"ModelType"] %in% "model control") & (out[,"InterventionType"] %in% "intervention control"), ]
        ctem = out[(out[,"ModelType"] %in% "model") & (out[,"InterventionType"] %in% "intervention control"), ]
        ttem = out[(out[,"ModelType"] %in% "model") & (out[,"InterventionType"] %in% "intervention"), ]
        
        # --------------------- SC -----------------
        if (type == "SC"){
          
          if (nrow(stem) == 0 | nrow(ctem) == 0) { next } # if only CT or ST, skip
          UniIntID = as.character(unique(ctem[,"InterventionID"]))
          
          for (k in 1:length(UniIntID)) { # within k-th Intervention type
            
            srec = stem[as.character(stem[,"InterventionID"]) %in% UniIntID[k],]     # nrow(srec) = 1
            if(nrow(srec) == 0) srec = stem[1,]
            
            crec = ctem[as.character(ctem[,"InterventionID"]) %in% UniIntID[k],]
            
            nr = nrow(crec)
            df = data.frame(matrix(nrow = nr, ncol = ncol(dat)+1), stringsAsFactors = F)
            colnames(df) = c("Author", "Year", "OutcomeId", "ch.s","ch.c", "ch.t", "Mean.S","Mean.C","Mean.T",
                             "SD.C", "SD.T", "ai", "bi", "ci", "di", "True.No.C", "No.T",
                             "ES","SE", colnames(dat[,(fix+1):ncol(dat), drop=FALSE]))
            
            df[,"Author"] = crec$Author
            df[,"Year"] = crec$Year
            df[,"OutcomeId"] = crec$OutcomeId
            df[,"ch.s"] = rep(NA, nr)
            df[,"ch.c"] = rep(srec$CohortId, nr)  # sham -> control
            df[,"ch.t"] = crec$CohortId           # control -> treatment
            df[,"True.No.C"] = rep(srec$NumberOfAnimals/nr, nr)
            df[,"No.T"] = crec$NumberOfAnimals
            df[,(fix+2):ncol(df)] = crec[,(fix+1):ncol(crec)]
            
            if (mea == "OR") {
              
              df[,"ai"] = crec$OutcomeResult # Number.Affectd.by.Outcome.Measure.in.Treatment.Group
              df[,"bi"] = crec$NumberOfAnimals - df$ai # Number.in.Treatment.Group - ai
              df[,"ci"] = srec$OutcomeResult # Number.Affectd.by.Outcome.Measure.in.Control.Group
              df[,"di"] = srec$NumberOfAnimals - df$ci # Number.in.Control.Group - ci
              df$ai[df$ai==0] = 0.5
              df$bi[df$bi==0] = 0.5
              df$ci[df$ci==0] = 0.5
              df$di[df$di==0] = 0.5
              
              df[,"ES"] = log(df[,"ai"]*df[,"di"]/(df[,"bi"]*df[,"ci"]))
              df[,"SE"] = 1/df[,"ai"] + 1/df[,"bi"] + 1/df[,"ci"] + 1/df[,"di"]
              
              df[,"Mean.S"] = df[,"Mean.C"] = df[,"Mean.T"]= rep(NA,nr)
              df[,"SD.C"] = df[,"SD.T"] = rep(NA,nr)
              
            } else { # mea == "NMD" or "SMD"
              
              df[,"Mean.S"] = rep(NA,nr)
              df[,"Mean.C"] = rep(srec$OutcomeResult, nr)
              df[,"Mean.T"]= crec$OutcomeResult
              
              if ( srec$ErrorType == "sd") {
                df[,"SD.C"] = rep(srec$OutcomeError, nr)
              } else { # srec$ErrorType == "sem"
                df[,"SD.C"] = rep(srec$OutcomeError*sqrt(srec$NumberOfAnimals), nr)
              } 
              
              for (r in 1:nr) { 
                if (crec$ErrorType[r] == "sd") {
                  df[r,"SD.T"] = crec$OutcomeError[r]
                } else { # crec$ErrorType[r] == "sem"
                  df[r,"SD.T"] = crec$OutcomeError[r]*sqrt(crec$NumberOfAnimals[r])
                }
              } 
              df[,"ai"] = df[,"bi"] = df[,"ci"] = df[,"di"] = rep(NA,nr)
              
              if (mea == "NMD") {
                df[,"ES"] = 100*(df[,"Mean.C"] - df[,"Mean.T"]) / df[,"Mean.C"]
                df[,"SE"] = (100*(df[,"SD.C"]/df[,"Mean.C"]))^2/df[,"True.No.C"] + (100*(df[,"SD.T"]/df[,"Mean.C"]))^2/df[,"No.T"]
                df[,"SE"] = sqrt(df[,"SE"])
              }
              if (mea == "SMD"){
                Sp = sqrt((df[,"True.No.C"]-1)*df[,"SD.C"]^2/(df[,"True.No.C"]+df[,"No.T"]-2) + (df[,"No.T"]-1)*df[,"SD.T"]^2/(df[,"True.No.C"]+df[,"No.T"]-2))
                df[,"ES"] = (df[,"Mean.C"]-df[,"Mean.T"])/Sp*(1-3/(4*(df[,"True.No.C"]+df[,"No.T"])-9))*crec[,"GreaterIsWorse"]
                df[,"SE"] = (df[,"True.No.C"]+df[,"No.T"])/(df[,"True.No.C"]*df[,"No.T"]) + df[,"ES"]^2/2/(df[,"True.No.C"]+df[,"No.T"]-3.94)
                df[,"SE"] = sqrt(df[,"SE"])
              }
            } # if (meth == "OR") 
            record = rbind(record, df) 
          } # for (k in 1:length(UniIntID)) 
        } # if (type == "SC") 
        
        # --------------------- SCT -----------------
        if (type == "SCT") {
          
          if (nrow(ttem) == 0 | nrow(ctem) == 0) { next } # if SC or ST, skip
          UniModID = as.character(unique(ttem$ModelID))
          
          for (k in 1:length(UniModID)) {
            
            crec = ctem[ctem$ModelID %in% UniModID[k],]   # nrow(crec) = 1
            trec = ttem[ttem$ModelID %in% UniModID[k],]
            nr = nrow(trec)
            
            df = data.frame(matrix(nrow = nr, ncol = ncol(dat)+1), stringsAsFactors=F)
            colnames(df) = c("Author", "Year", "OutcomeId", "ch.s","ch.c", "ch.t", "Mean.S","Mean.C","Mean.T",
                             "SD.C", "SD.T", "ai", "bi", "ci", "di", "True.No.C", "No.T",
                             "ES","SE", colnames(dat[,(fix+1):ncol(dat), drop=FALSE]))
            
            df[,"Author"] = trec$Author
            df[,"Year"] = trec$Year
            df[,"OutcomeId"] = crec$OutcomeId
            if (nrow(stem) == 0) {df[,"ch.s"] = rep(NA, nr)}
            if (nrow(stem) == 1) {df[,"ch.s"] = rep(stem$CohortId, nr)}
            df[,"ch.c"] = rep(crec$CohortId, nr)  # sham -> control
            df[,"ch.t"] = trec$CohortId           # control -> treatment
            df[,"True.No.C"] = rep(crec$NumberOfAnimals/nr, nr)
            df[,"No.T"] = trec$NumberOfAnimals
            df[,(fix+2):ncol(df)] = trec[,(fix+1):ncol(trec)]
            
            if (mea == "OR") {
              
              df[,"ai"] = trec$OutcomeResult # Number.Affectd.by.Outcome.Measure.in.Treatment.Group
              df[,"bi"] = trec$NumberOfAnimals - df$ai # Number.in.Treatment.Group - ai
              df[,"ci"] = crec$OutcomeResult # Number.Affectd.by.Outcome.Measure.in.Control.Group
              df[,"di"] = crec$NumberOfAnimals - df$ci # Number.in.Control.Group - ci
              df$ai[df$ai==0] = 0.5
              df$bi[df$bi==0] = 0.5
              df$ci[df$ci==0] = 0.5
              df$di[df$di==0] = 0.5
              
              df[,"ES"] = log(df[,"ai"]*df[,"di"]/(df[,"bi"]*df[,"ci"]))
              df[,"SE"] = 1/df[,"ai"] + 1/df[,"bi"] + 1/df[,"ci"] + 1/df[,"di"]
              
              df[,"Mean.S"] = df[,"Mean.C"] = df[,"Mean.T"]= rep(NA,nr)
              df[,"SD.C"] = df[,"SD.T"] = rep(NA,nr)
              
            } else { # mea == "NMD" or "SMD"
              
              if (nrow(stem) == 0) {
                df[,"Mean.S"] = NA
                if (is.na(trec[,"InferredSham"])){
                  df[,"Mean.S"] = rep(0, nr)
                } else {
                  df[,"Mean.S"] = trec[,"InferredSham"]
                }
              }
              if (nrow(stem) == 1) {
                df[,"Mean.S"] = rep(stem$OutcomeResult, nr)
                df[is.na(df[,"Mean.S"]), "Mean.S"] = 0
              } 
              
              df[,"Mean.C"] = rep(crec$OutcomeResult, nr)
              df[,"Mean.T"]= trec$OutcomeResult
              
              if (crec$ErrorType == "sd") {
                df[,"SD.C"] = rep(crec$OutcomeError, nr)
              } else { # crec$ErrorType == "sem"
                df[,"SD.C"] = rep(crec$OutcomeError*sqrt(crec$NumberOfAnimals), nr)
              } 
              for (r in 1:nr) { 
                if (trec$ErrorType[r] == "sd") {
                  df[r,"SD.T"] = trec$OutcomeError[r]
                } else { # trec$ErrorType[r] == "sem"
                  df[r,"SD.T"] = trec$OutcomeError[r]*sqrt(trec$NumberOfAnimals[r])
                }
              } 
              
              df[,"ai"] = df[,"bi"] = df[,"ci"] = df[,"di"] = rep(NA,nr)
              
              if (mea == "NMD") {
                df[,"ES"] = 100*(df[,"Mean.C"] - df[,"Mean.T"]) / (df[,"Mean.C"] - df[,"Mean.S"])
                df[,"SE"] = (100*(df[,"SD.C"]/(df[,"Mean.C"]-df[,"Mean.S"])))^2/df[,"True.No.C"] + (100*(df[,"SD.T"]/(df[,"Mean.C"]-df[,"Mean.S"])))^2/df[,"No.T"]
                df[,"SE"] = sqrt(df[,"SE"])
              }
              if (mea == "SMD"){
                Sp = sqrt((df[,"True.No.C"]-1)*df[,"SD.C"]^2/(df[,"True.No.C"]+df[,"No.T"]-2) + (df[,"No.T"]-1)*df[,"SD.T"]^2/(df[,"True.No.C"]+df[,"No.T"]-2))
                df[,"ES"] = (df[,"Mean.C"]-df[,"Mean.T"])/Sp*(1-3/(4*(df[,"True.No.C"]+df[,"No.T"])-9))*trec[,"GreaterIsWorse"]
                df[,"SE"] = (df[,"True.No.C"]+df[,"No.T"])/(df[,"True.No.C"]*df[,"No.T"]) + df[,"ES"]^2/2/(df[,"True.No.C"]+df[,"No.T"]-3.94)
                df[,"SE"] = sqrt(df[,"SE"])
              }
            } # if (meth == "OR") 
            record = rbind(record, df) 
          } # for (k in 1:length(UniModID))
        } # if (type == "SCT")
      } # for (j in 1:length(OutList))
    } # for (i in 1:length(ExpList))
    
    return(record)
  })
  
  
  # ----------------- Nest variables ---------------------------------
  NestVar <- reactive({
    
    pre = record()
    pre[, "ch.s"] = paste0(pre[,"ch.s"], pre[,"ch.c"], pre[,"ch.t"])
    pre = pre[, !(colnames(pre) %in% c("ch.c", "ch.t"))]
    
    flag = matrix(0, nrow = 1, ncol = ncol(pre)-17)
    colnames(flag) = colnames(pre[,18:ncol(pre)])
    flag = data.frame(flag)
    
    corlist = split(pre, pre[,"ch.s"], drop = T)
    for (i in 1:length(corlist)){
      temp = corlist[[i]]
      nc = ncol(temp)-17
      for (j in 1:nc){
        if (length(unique(temp[,j+17])) > 1) {
          flag[colnames(temp)[j+17]] = 1
        }
      }
    }
    
    NestVar = colnames(flag)[flag==1]
    return(NestVar)
    
  })
  
  # ----------------- Nest way ---------------------------------
  output$NestWay <- renderUI({
    
    nest_name = NestVar()
    fluidRow(
      lapply(1:length(nest_name), function(i) {
        column(3,
               selectInput(inputId = paste0("nest_",i), label = nest_name[i],
                           choices = list("First" = "1st",
                                          "Last" = "last",
                                          "Count" = "count",
                                          "Sum" = "sum",
                                          "Average" = "mean",
                                          "Minimum" = "min",
                                          "Maximum" = "max",
                                          "Variance" = "var",
                                          "Standard deviation" = "sd",
                                          "Standard error" = "se"
                           ),
                           multiple = F, selected = "1st")
        ) # column
      }) # lapply(1:length(nest_name), function(i)
    ) # fluidRow
    
  })
  
  # ----------------- Nest ---------------------------------
  nest <- reactive({
    
    pre = record()
    pre[, "ch.s"] = paste0(pre[,"ch.s"], pre[,"ch.c"], pre[,"ch.t"])
    pre = pre[, !(colnames(pre) %in% c("ch.c", "ch.t"))]
    
    nest = c()
    corlist <- split(pre, pre[,"ch.s"], drop = T)
    
    nest_name = NestVar()
    all_name = colnames(pre[,20:ncol(pre)])
    rest_name = all_name[!(all_name %in% nest_name)]
    
    for (i in 1:length(corlist)) {
      temp = corlist[[i]] # i-th cohort group
      UniOutID = as.character(unique(temp[,"OutcomeId"]))
      for (j in 1:length(UniOutID)){
        out = temp[temp[,"OutcomeId"] %in% UniOutID[j],] # in i-th Cohort group with same j-th OutcomeId
        df = data.frame(matrix(nrow = 1, ncol = 6+ncol(pre)-17), stringsAsFactors = F)
        colnames(df) = c("Author", "Year", "True.No.C", "No.T", "ES","SE",
                         colnames(pre[,18:ncol(pre)]))
        
        df[1,"Author"] = as.character(temp[1,"Author"])
        df[1,"Year"] = temp[1,"Year"]
        df[1,"True.No.C"] = temp[1,"True.No.C"]
        df[1,"No.T"] = temp[1,"No.T"]
        
        yi = temp[,"ES"]
        wi = 1/(temp[,"SE"]^2)
        df[,"ES"] = sum(yi*wi)/sum(wi)
        df[,"SE"] = sqrt(nrow(temp)/sum(wi))
        
        
        for (r in 1:length(rest_name)) {
          df[,rest_name[r]] = temp[1,rest_name[r]]
        }
        
        for (k in 1:length(nest_name)) {
          
          nest_way = input[[paste0("nest_", k)]]
          
          if (nest_way == "1st") {
            df[,nest_name[k]] = temp[1,nest_name[k]]
          }
          if (nest_way == "last") {
            df[,nest_name[k]] = temp[nrow(temp),nest_name[k]]
          }
          if (nest_way == "count") {
            df[,nest_name[k]] = nrow(temp)
          }
          if (nest_way == "sum") {
            df[,nest_name[k]] = sum(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "mean") {
            df[,nest_name[k]] = mean(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "min") {
            df[,nest_name[k]] = min(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "max") {
            df[,nest_name[k]] = max(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "var") {
            df[,nest_name[k]] = var(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "sd") {
            df[,nest_name[k]] = sd(temp[,nest_name[k]], na.rm = T)
          }
          if (nest_way == "se") {
            df[,nest_name[k]] = sd(temp[,nest_name[k]], na.rm = T)/sqrt(nrow(temp))
          }
        } # for (k in 1:length(nest_name))
        
        nest = rbind(nest, df)
        
      } # for (j in 1:length(UniOutID))
    } # for (i in 1:length(corlist))
    
    return(nest)
    
  })
  
  
  # ----------------- DT ----------------------
  output$DT <- renderDataTable({
    
    if (is.null(dataset())) { stop("No data submitted.") }
    
    if (input$DataType == "pre") {
      dat <- record()
      dat = dat [, !(colnames(dat) %in% c("ch.s", "ch.c", "ch.t", "OutcomeId"))]
      ndigit = 2
      
      dat[,"SD.C"] = round(dat[,"SD.C"], ndigit)
      dat[,"SD.T"] = round(dat[,"SD.T"], ndigit)
      dat[,"ai"] = round(dat[,"ai"], ndigit)
      dat[,"bi"] = round(dat[,"bi"], ndigit)
      dat[,"ci"] = round(dat[,"ci"], ndigit)
      dat[,"di"] = round(dat[,"di"], ndigit)
      dat[,"True.No.C"] = round(dat[,"True.No.C"], ndigit)
      dat[,"No.T"] = round(dat[,"No.T"], ndigit)
      dat[,"ES"] = round(dat[,"ES"], ndigit)
      dat[,"SE"] = round(dat[,"SE"], ndigit)
      
      DT = dat
    }
    
    if (input$DataType == "nest") {
      dat = nest()
      ndigit = 2
      
      dat[,"True.No.C"] = round(dat[,"True.No.C"], ndigit)
      dat[,"No.T"] = round(dat[,"No.T"], ndigit)
      dat[,"ES"] = round(dat[,"ES"], ndigit)
      dat[,"SE"] = round(dat[,"SE"], ndigit)
      
      DT = dat
    }
    
    return (DT)
  }, options = list(pageLength = 10))   
  
  
  output$DownTable <- downloadHandler(
    filename = function() {
      paste("Table", input$FileType, sep = ".")
    },
    content = function(file) {
      if (input$DataType == "pre"){ 
        dat <- record()
        dat = dat[, !(colnames(dat) %in% c("ch.s", "ch.c", "ch.t", "OutcomeId"))]
        DT <- dat 
      }
      if (input$DataType == "nest"){ 
        DT <- nest() 
      }
      sep <- switch(input$FileType, "csv" = ",", "tsv" = "\t")
      write.table(DT, file, sep = sep, row.names = F)
    }
  )
  
  
  # ------------------- Tell SyRF continuous variables --------------------
  output$CheckVar <- renderUI({
    
    if (input$DataType == "pre") { 
      dat <- record()
      use <- names(dat)[20:ncol(dat)]
    }
    if (input$DataType == "nest") {
      dat <- nest()
      use <- names(dat)[7:ncol(dat)]
    }
    
    selectInput(inputId = "check_var", label = "Tell SyRF which variables are continuous",
                choices = as.list(use), multiple = T, selectize = F) 
    
  })
  
  # ----------------- Meta-analysis ----------------------
  output$GlobalOutput <- renderPrint({
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest() }
    if (is.null(yivi)) { stop("No data submitted or select wrong model.") }
    
    # summary(metagen(Effect.Size, Standard.Error, data=yv, method.tau=input$meth))
    # summary(metagen(TE = ES, seTE = SE, data=yivi, method.tau="REML"))
    summary(metagen(TE = ES, seTE = SE, data = yivi,
                    method.tau = input$HetEstimator))
  })
  
  # ----------------- Forest plot ----------------------
  Fplot <- function(){
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest()}
    # if (is.null(yivi)) { stop("No data submitted.") }
    
    res <- metagen(TE = ES, seTE = SE, data = yivi,
                   method.tau = input$HetEstimator, comb.fixed = F)
    
    if (input$ForOrder == "") { sortlab = NULL }
    if (input$ForOrder == "ine") { sortlab = res$TE }
    if (input$ForOrder == "inw") { sortlab = res$w.random }
    if (input$ForOrder == "iny") { sortlab = yivi$Year }
    if (input$ForOrder == "dee") { sortlab = -res$TE }
    if (input$ForOrder == "dew") { sortlab = -res$w.random }
    if (input$ForOrder == "dey") { sortlab = -yivi$Year }
    
    
    if (input$ShowWeight == T) {
      right_cols = c("effect", "ci", "w.random")
      right_labs = c("Effect","95% CI", "Weight")
    } else {
      right_cols = c("effect", "ci")
      right_labs = c("Effect","95% CI")
    }
    
    forest(res, hetstat = F,
           studlab = paste0(yivi[,"Author"], ", ", yivi[,"Year"]),
           leftcols = c("studlab"),
           rightcols = right_cols, rightlabs = right_labs,
           sortvar = sortlab,# Year, w.random, TE
           
           xlab = input$ForXlab,
           col.square = input$ForSqCol,
           col.diamond = input$ForDiaCol,
           plotwidth = paste0(input$ForWidth, "cm"),
           colgap.forest.left = paste0(input$GapLeft, "cm"),
           colgap.forest.right = paste0(input$GapRight, "cm"),
           col.square.lines = "black",
           col.diamond.lines = "black",
           col.i="black"
    )
  }
  
  FWinW <- reactive({ input$ForWinWidth })
  FWinH <- reactive({ input$ForWinHeight })
  
  output$ForestPlot <- renderPlot({
    Fplot()
  }, width = FWinW, height = FWinH)
  
  output$DownForest <- downloadHandler(
    filename = function() {
      paste("forest-", Sys.time(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = FWinW(), height = FWinH()) # open the png device
      Fplot()
      dev.off()  # turn the device off
    }
  )
  
  # ----------------- Select variable for Het analysis ----------------------
  # for stratified (one discrete variable)
  output$SubVar <- renderUI({
    
    if (input$DataType == "pre") {
      dat <- record()
      use <- names(dat)[20:ncol(dat)]
    }
    if (input$DataType == "nest") {
      dat <- nest()
      use <- names(dat)[7:ncol(dat)]
    }
    
    dname = use[!use %in% input$check_var]  # Name of all discrete variables
    
    selectInput(inputId = "subvar",
                label = "Select discrete variable for stratified meta-analysis",
                choices = as.list(dname), multiple = F, selectize = F)
  })
  
  # for meta-regression (multi continuous variables)
  output$RegContVar <- renderUI({
    
    cname = input$check_var                 # Name of all continuous varoables
    selectInput(inputId = "cbox",
                label = "Select continuous variables for meta-regression",
                choices = as.list(cname), multiple = T, selectize = F)
  })
  
  # for meta-regression (multi discrete variables)
  output$RegDiscVar <- renderUI({
    
    if (input$DataType == "pre") {
      dat <- record()
      use <- names(dat)[20:ncol(dat)]
    }
    if (input$DataType == "nest") {
      dat <- nest()
      use <- names(dat)[7:ncol(dat)]
    }
    
    dname = use[!use %in% input$check_var]  # Name of all discrete variables
    selectInput(inputId = "dbox",
                label = "Select discrete variables for meta-regression",
                choices = as.list(dname), multiple = T, selectize = F)
  })
  
  # ---------------------- Het (stratified) ------------------------
  output$SubOutput <- renderPrint({
    
    # if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest() }
    
    
    if (input$DataType == "pre") { 
      # stop("Pre-nested data shouldn't be used for heterogeneity analysis. 
      #      Please select nested data in the first tab.")
      yivi <- record()
    }
    
    summary(metagen(TE = ES, seTE = SE, data = yivi,
                    method.tau = input$HetEstimator,
                    byvar = yivi[,input$subvar]))
  })
  
  # ---------------------- Het (meta-regression) ------------------------
  output$RegOutput <- renderPrint({
    
    # if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest()}
    
    if (input$DataType == "pre") { 
      # stop("Pre-nested data shouldn't be used for heterogeneity analysis. 
      #      Please select nested data in the first tab.")
      yivi <- record()
    }
    
    dlab = input$dbox   # discrete variables selected
    clab = input$cbox   # continuous variables selected
    
    if(length(clab) == 0 & length(dlab) == 0) { stop("No variable selected") }
    if(length(clab) == 0 & length(dlab) != 0) { # only discrete variables
      mo = paste("~",paste(paste("factor(", paste(dlab, collapse=")+factor(")),")"), collapse="")
    }
    
    if(length(clab) != 0 & length(dlab) == 0){ # only continuous variables
      mo = paste("~",paste(clab, collapse="+"), collapse="")
    }
    
    if(length(clab) != 0 & length(dlab) != 0){ # discrete and continuous
      form.d = paste(paste("+factor(", paste(dlab, collapse=")+factor(")),")")
      form.c = paste("~",paste(clab, collapse="+"), collapse="")
      mo = paste(form.c, form.d, collapse="+")
    }
    
    rma(yi = ES, vi = (SE)^2, mods = as.formula(mo),
        data = yivi, method = input$HetEstimator)
  })
  
  
  # ------------------- Select box for Het plot --------------------
  # for sub forest plot (one discrete variable)
  output$SubPlotVar <- renderUI({
    
    dat <- record()
    allname <- names(dat)[20:ncol(dat)]
    dname = allname[!allname %in% input$check_var]
    
    selectInput(inputId = "dis.forest",
                label = "Select a variable for forest plot",
                choices = as.list(dname), selectize = F)
  })
  
  # for meta-regression plot
  output$RegPlotVar <- renderUI({
    
    cname = input$check_var
    selectInput(inputId = "con.regplot",
                label = "Select a variable for meta-regression plot",
                choices = as.list(cname), selectize = F)
  })
  
  
  # ---------------------- Het plot (stratified: subforest plot) ------------------------
  Fplot2 <- function(){
    
    if (length(input$dis.forest) == 0) { stop("No discrete variable selected.") }
    
    # if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest() }
    
    if (input$DataType == "pre") { 
      # stop("Pre-nested data shouldn't be used for heterogeneity analysis. 
      #      Please select nested data in the first tab.")
      yivi <- record()
    }
    
    res <- metagen(TE = ES, seTE = SE, data = yivi,
                   method.tau = input$HetEstimator,
                   byvar = yivi[,input$dis.forest])
    
    if (input$ShowFixWeight == T & input$ShowRandWeight == T) {
      right_cols = c("effect", "ci", "w.fixed", "w.random")
      right_labs = c("Effect", "95% CI", "W-fixed", "W-random")
    }
    if (input$ShowFixWeight == T & input$ShowRandWeight == F) {
      right_cols = c("effect", "ci", "w.fixed")
      right_labs = c("Effect", "95% CI", "W-fixed")
    }
    if (input$ShowFixWeight == F & input$ShowRandWeight == T) {
      right_cols = c("effect", "ci", "w.random")
      right_labs = c("Effect", "95% CI", "W-random")
    }
    if (input$ShowFixWeight == F & input$ShowRandWeight == F) {
      right_cols = c("effect", "ci")
      right_labs = c("Effect", "95% CI")
    }
    
    forest(res, hetstat = F,
           bylab = input$dis.forest,
           print.subgroup.labels = T,
           studlab = paste0(yivi$Author, ", ", yivi$Year),
           leftcols = c("studlab"),
           rightcols = right_cols,
           rightlabs = right_labs,
           
           xlab = input$SubXlab,
           col.square = input$SubSqCol,
           col.diamond = input$SubDiaCol,
           plotwidth = paste0(input$SubWidth, "cm"),
           colgap.forest.left = paste0(input$SubGapLeft, "cm"),
           colgap.forest.right = paste0(input$SubGapRight, "cm"),
           
           col.square.lines = "black",
           col.diamond.lines = "black",
           col.i = "black"
    )
  }
  
  
  SWinW <- reactive({ input$SubWinWidth })
  SWinH <- reactive({ input$SubWinHeight })
  
  output$SubForest <- renderPlot({
    Fplot2()
  }, width = SWinW, height = SWinH)
  
  output$DownSubForest <- downloadHandler(
    filename = function() {
      paste("forest-", Sys.time(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = SWinW(), height = SWinH()) # open the png device
      Fplot2()
      dev.off()  # turn the device off
    }
  )
  
  # ----------------- Het plot (meta-regression) ----------------------
  output$RegPlot <- renderPlotly({
    
    if (length(input$con.regplot)==0 | length(input$check_var)==0) { stop("No continuous variable selected.") }
    
    # if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest()}
    if (input$DataType == "pre") { 
      # stop("Pre-nested data shouldn't be used for heterogeneity analysis. 
      #      Please select nested data in the first tab.")
      yivi <- record()
    }
    
    lab = input$con.regplot          # Name of the continuous variable
    Cont = yivi[,lab]
    
    res <- rma(yi = ES, vi = SE^2, mods = ~Cont, data = yivi,
               method = input$HetEstimator)
    
    Size = 1/(yivi[,"SE"])^2
    yfit = Cont*res$b[2]+res$b[1]
    
    if (input$RegXlab == "type the xlabel") {
      reg_xlab = lab
    } else {
      reg_xlab = input$RegXlab
    }
    reg_ylab = input$RegYlab
    
    plot_ly() %>%
      add_trace(data = yivi,
                x = ~Cont,
                y = ~ES,
                size = ~Size,
                mode = "markers",
                showlegend = F,
                name = "Observed",
                
                hoverinfo = 'text',
                text = ~paste('Effect: ', round(ES,2),
                              '</br>SE: ', round(yivi$SE,2),
                              '</br>Author: ', Author,
                              '</br>Year: ', Year
                ),
                marker = list(color = input$RegPtCol)
      ) %>%
      
      add_trace(data = yivi,
                x = ~Cont,
                y = yfit,
                mode = "lines",
                showlegend = F,
                name = "Meta-regression",
                line = list(color = input$RegLCol)
      ) %>%
      layout(
        xaxis = list(title = reg_xlab,
                     showgrid = F,
                     showline = T,
                     zeroline = T,
                     zerolinecolor = toRGB("gray90"),
                     zerolinewidth = 1),
        
        yaxis = list(title = reg_ylab,
                     showline = T,
                     zeroline = T,
                     zerolinecolor = toRGB("gray90"),
                     zerolinewidth = 1),
        autosize = F,
        width = input$RegWidth,
        height = input$RegHeight
      ) # layout
  })
  
  
  # ------------------- Select variables for bar plot --------------------
  output$BarVar <- renderUI({
    
    if (input$DataType == "pre") {
      dat <- record()
      allname <- names(dat)[20:ncol(dat)]
    }
    if (input$DataType == "nest") {
      dat <- nest()
      allname <- names(dat)[7:ncol(dat)]
    }
    
    dname = allname[!allname %in% input$check_var]
    selectInput(inputId = "dis.bar",
                label = "Select variables for bar plot",
                choices = as.list(dname), multiple = F, selectize = F)
  })
  
  
  # ----------------- Bar plot ----------------------
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      for (i in 1:numPlots) {
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  Bplot <- function(){
    
    lab = input$dis.bar
    if (is.null(lab)) { stop("No variable selected") }
    nbar = length(lab)
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest() }
    
    yivi$ni = yivi[,"True.No.C"] + yivi[,"No.T"]
    
    # stratified bar plot
    if (input$HetMethod == "sub"){
      
      plots <- list()  # new plots list
      for (j in 1:nbar){
        
        mod <- metagen(TE = ES, seTE = SE, data = yivi, comb.fixed = F,
                       method.tau = input$HetEstimator,
                       byvar = yivi[,lab[[j]]])
        
        tem = aggregate(yivi$ni, by = list(gro=yivi[,lab[[j]]]), FUN = sum)
        temp = tem
        tlab = as.matrix(mod$bylevs)
        
        for (t in 1:nrow(temp)) {
          temp[t,] = tem[tem[,1] == tlab[t],]
        }
        
        df = data.frame(xlabel = temp[,1],
                        effect = mod$TE.random.w,
                        w = sqrt(temp[,2])/sum(sqrt(temp[,2])),
                        c_low = mod$lower.random.w,
                        c_up = mod$upper.random.w)
        
        bar_ylab = input$BarYlab
        if (input$BarTitle == "type the title") {
          bar_title = lab[[j]]
        } else {
          bar_title = input$BarTitle
        }
        
        if (input$BarYmin == "") {
          y_min = min(0, df$c_low, mod$lower.random)
        } else {
          y_min = as.numeric(input$BarYmin)
        }
        if (input$BarYmax == "") {
          y_max = max(0, df$c_up, mod$upper.random)
        } else {
          y_max = as.numeric(input$BarYmax)
        }
        
        p = ggplot(df, aes(x = xlabel,
                           y = effect,
                           width = w)) +
          
          geom_col(fill = "white", colour = "black") +
          geom_rect(aes(xmin = -Inf,
                        xmax = Inf,
                        ymin = mod$lower.random,
                        ymax = mod$upper.random),
                    fill = "grey90", alpha = 0.2) +
          
          geom_linerange(aes(ymin = c_low,
                             ymax = c_up),
                         data = df) +
          
          geom_hline(yintercept = 0, linetype = 2) +
          scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) +
          
          
          labs(x = "",y = bar_ylab) +
          ggtitle(bar_title) +
          theme(axis.line = element_line(colour = "black"),
                axis.text = element_text(size = 12),
                plot.title = element_text(size = input$BarTitleSize,
                                          hjust = 0.5,
                                          vjust = 0.5),
                
                axis.title = element_text(size = input$BarYlabSize,
                                          hjust = 0.5,
                                          vjust = 0.5),
                
                axis.text.x = element_text(size = input$BarLabSize,
                                           hjust = input$BarLabPos,
                                           angle = input$BarLabAngle),
                
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank()
          )
        
        plots[[j]] = p
      }  # for (j in 1:nbar)
      multiplot(plotlist = plots, cols = 2)
    }
    
    
    # meta-regression bar plot
    if (input$HetMethod == "reg"){
      
      plots <- list()  # new plots list
      for (j in 1:nbar){
        
        # Split into subgroups
        mylist <- split(yivi, yivi[,lab[[j]]], drop = T)
        
        
        # Do meta-regression
        res <- rma(yi = ES, vi = SE^2, mods = ~factor(yivi[,lab[[j]]])-1,
                   data = yivi, method = input$HetEstimator)
        glob <- rma(yi = ES, vi = SE^2, data = yivi,
                    method = input$HetEstimator)
        
        lev = rownames(res$b) # names for all the levels of one variable
        tem = aggregate(yivi$ni, by = list(gro=yivi[,lab[[j]]]), FUN = sum)
        
        df = data.frame(xlabel = substr(lev, 25, nchar(lev)),
                        y = res$b,
                        low = res$ci.lb,
                        up = res$ci.ub,
                        n = sqrt(tem[,2]))
        
        
        bar_ylab = input$BarYlab
        if (input$BarTitle == "type the title") {
          bar_title = lab[[j]]
        } else {
          bar_title = input$BarTitle
        }
        
        if (input$BarYmin == "") {
          y_min = min(0, df$low, glob$ci.lb)
        } else {
          y_min = as.numeric(input$BarYmin)
        }
        
        if (input$BarYmax == "") {
          y_max = max(0, df$up, glob$ci.ub)
        } else {
          y_max = as.numeric(input$BarYmax)
        }
        
        
        p = ggplot(df, aes(x = xlabel,
                           y = y,
                           width = n/sum(n))) +
          
          geom_col(fill = "white", colour = "black") +
          
          geom_rect(aes(xmin = -Inf,
                        xmax = Inf,
                        ymin = glob$ci.lb,
                        ymax = glob$ci.ub),
                    fill = "grey90", alpha = 0.2) +
          
          geom_linerange(aes(ymin = low,
                             ymax = up),
                         data = df) +
          
          geom_hline(yintercept = 0, linetype = 2) +
          
          scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) +
          
          labs(x = "", y = bar_ylab) +
          ggtitle(bar_title) +
          theme(axis.line = element_line(colour = "black"),
                
                axis.text = element_text(size = 12),
                
                plot.title = element_text(size = input$BarTitleSize,
                                          hjust = 0.5,
                                          vjust = 0.5),
                
                axis.title = element_text(size = input$BarYlabSize,
                                          hjust = 0.5,
                                          vjust = 0.5), # y
                
                axis.text.x = element_text(size = input$BarLabSize,
                                           hjust = input$BarLabPos,
                                           angle = input$BarLabAngle),  # bar label
                
                
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank()
          )
        plots[[j]] = p
      }  # for (j in 1:nbar)
      multiplot(plotlist = plots, cols = 2)
    } # if (input$het == "reg")
  }
  
  
  BarW <- reactive({ input$BarWidth })
  BarH <- reactive({ input$BarHeight })
  
  output$BarPlot <- renderPlot({
    Bplot()
  }, width = BarW, height = BarH)
  
  output$DownBar <- downloadHandler(
    filename = function() {
      paste("Bar-", Sys.time(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = BarW(), height = BarH()) # open the png device
      Bplot()
      dev.off()  # turn the device off
    }
  )
  
  
  
  Binfo <- function(){
    
    lab = input$dis.bar
    if (is.null(lab)) { stop("No variable selected") }
    nbar = length(lab)
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { yivi <- nest() }
    
    yivi$ni = yivi[,"True.No.C"] + yivi[,"No.T"]
    
    # stratified bar info
    if (input$HetMethod == "sub"){
      
      mod <- metagen(TE = ES, seTE = SE, data = yivi, comb.fixed = F,
                     method.tau = input$HetEstimator,
                     byvar = yivi[,lab])
      
      tem = aggregate(yivi$ni, by = list(gro=yivi[,lab]), FUN = sum)
      temp = tem
      tlab = as.matrix(mod$bylevs)
      
      for (t in 1:nrow(temp)) {
        temp[t,] = tem[tem[,1] == tlab[t],]
      }
      
      ndig = 2
      
      df = data.frame(Group = temp[,1],
                      ES = round(mod$TE.random.w, ndig),
                      ci.low = round(mod$lower.random.w, ndig),
                      ci.up = round(mod$upper.random.w, ndig),
                      No.Animals = temp[,2],
                      No.Studies = mod$k.w,
                      stringsAsFactors = F)
      df = rbind(df,
                 c("Global", 
                   round(mod$TE.random, ndig),
                   round(mod$lower.random, ndig),
                   round(mod$upper.random, ndig), 
                   sum(yivi[,"ni"]), mod$k)
      )
    } # if (input$HetMethod == "sub")
    
    # meta-regression bar info
    if (input$HetMethod == "reg"){
      
      # Split into subgroups
      mylist <- split(yivi, yivi[,lab], drop = T)
      # Do meta-regression
      res <- rma(yi = ES, vi = SE^2, mods = ~factor(yivi[,lab])-1,
                 data = yivi, method = input$HetEstimator)
      glob <- rma(yi = ES, vi = SE^2, data = yivi,
                  method = input$HetEstimator)
      
      lev = rownames(res$b) # names for all the levels of one variable
      
      ndig = 2
      df = data.frame(Group = substr(lev, 20, nchar(lev)),
                      ES = round(res$b, ndig),
                      ci.low = round(res$ci.lb ,ndig),
                      ci.up = round(res$ci.ub, ndig),
                      stringsAsFactors = F)
      df = rbind(df,
                 c("Global", 
                   round(glob$b, ndig), 
                   round(glob$ci.lb, ndig), 
                   round(glob$ci.ub, ndig))
      )
    } # if (input$het == "reg")
    return(df)
  } 
  
  output$BarInfoTable <- renderTable({
    Binfo()
  })
  
  # ---------------- Trim-and-Fill -------------------------
  output$TafOutput <- renderPrint({
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { 
      stop("Nested data shouldn't be used for trim and fill analysis. 
           Please select pre-nested data in the first tab.")
      # yivi <- nest() 
    }
    # if (is.null(yivi)) { stop("No data submitted.") }
    
    res <- rma(yi = ES, vi = SE^2, ni = True.No.C+No.T,
               data = yivi, method = input$HetEstimator)
    trimfill(res, estimator = "L0", side = input$TafSide)
  })
  
  
  # ----------------- Funnel plot ----------------------
  Funplot <- function(){
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { stop() }
    # if (is.null(yivi)) { stop("No data submitted.") }
    
    
    res <- rma(yi = ES, vi = SE^2, ni = True.No.C+No.T,
               data = yivi, method = input$HetEstimator)
    taf <- trimfill(res, estimator="L0", side = input$TafSide)
    
    xmin = min(min(taf$yi)-0.1*(max(taf$yi)-min(taf$yi)), 0)
    xmax = max(max(taf$yi)+0.1*(max(taf$yi)-min(taf$yi)), 0)
    
    
    if (input$TafYaxis == "seinv"){
      ymax = max(max(1/sqrt(taf$vi))+0.1*(max(1/sqrt(taf$vi))-min(1/sqrt(taf$vi))), 0)
      if (input$TafFill == "Yes"){
        funnel(taf,
               xlim = c(xmin, xmax),
               ylim = c(0.000000001, ymax),
               back = "white",
               refline = res$b,
               yaxis = input$TafYaxis,
               xlab = input$TafXlab,
               ylab = input$TafYlab)
      } else {
        funnel(res,
               xlim = c(xmin, xmax),
               ylim = c(0.000000001, ymax),
               back = "white",
               yaxis = input$TafYaxis,
               xlab = input$TafXlab,
               ylab = input$TafYlab)
      }
    }
    
    if (input$TafYaxis == "sqrtninv"){
      ymax = max(max(1/sqrt(taf$ni))+0.1*(max(1/sqrt(taf$ni))-min(1/sqrt(taf$ni))), 0)
      
      if (input$TafFill == "Yes"){
        funnel(taf,
               xlim = c(xmin, xmax),
               ylim = c(0.000000001, ymax),
               back = "white",
               refline = res$b,
               yaxis = input$TafYaxis,
               xlab = input$TafXlab,
               ylab = input$TafYlab)
      } else {
        funnel(res,
               xlim = c(xmin, xmax),
               ylim = c(0.000000001, ymax),
               back = "white",
               yaxis = input$TafYaxis,
               xlab = input$TafXlab,
               ylab = input$TafYlab)
      }
    }
    abline(v = taf$b, lty = 2)
  }
  
  FnlW <- reactive({ input$FunnelWidth })
  FnlH <- reactive({ input$FunnelHeight })
  
  
  output$FunnelPlot <- renderPlot({
    Funplot()
  },  width = FnlW, height = FnlH)
  
  output$DownFunnel <- downloadHandler(
    filename = function() {
      paste("funnel-", Sys.time(), ".png", sep="")
    },
    content = function(file) {
      png(file, width = FnlW(), height = FnlH()) # open the png device
      Funplot()
      dev.off()  # turn the device off
    }
  )
  
  # ---------------- Egger's regression -------------------------
  output$EggerOutput <- renderPrint({
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { stop() }
    # if (is.null(yivi)) { stop("No data submitted.") }
    
    Standardised.effect <- yivi[,"ES"]/yivi[,"SE"]
    Precision <- 1/yivi[,"SE"]
    egger <- lm(Standardised.effect ~ Precision)
    summary(egger)
  })
  
  
  # ----------------- Egger's regression plot ----------------------
  output$EggerPlot <- renderPlotly({
    
    if (input$DataType == "pre") { yivi <- record() }
    if (input$DataType == "nest") { 
      stop("Nested data shouldn't be used for trim and fill analysis. 
           Please select pre-nested data in the first tab.")
      # yivi <- nest() 
    }
    # if (is.null(yivi)) { stop("No data submitted.") }
    
    Standardised.effect <- yivi[,"ES"]/yivi[,"SE"]
    Precision <- 1/yivi[,"SE"]
    egger <- lm(Standardised.effect ~ Precision)
    
    xnew = seq(0, to = max(Precision)+(max(Precision)-min(Precision))*0.1,
               len = length(Precision))
    ynew = egger$coefficients[1] + xnew*egger$coefficients[2]
    pnew <- predict(egger, newdata = data.frame(Precision=xnew),
                    interval = "confidence", level = 0.95)
    
    
    # RibCol = col2rgb(input$EggerRibCol)
    
    plot_ly() %>%
      
      add_ribbons(x = ~xnew,
                  ymin = pnew[,"lwr"],
                  ymax = pnew[,"upr"],
                  line = list(color = "transparent"),
                  fillcolor = col2rgb(input$EggerRibCol),
                  opacity = 0.2,
                  showlegend = F, name = "95% confidence") %>%
      
      add_trace(x = ~xnew,
                y = ynew,
                mode = "lines",
                line = list(color = input$EggerLCol, width = 1.5),
                showlegend = F, name = "Egger's fit") %>%
      
      # add_trace(x = ~xnew,
      #           y = pnew[,"lwr"],
      #           mode = "lines",
      #           line = list(color = input$EggerRibCol, width = 1.5),
      #           showlegend = F, name = "95% lower CI") %>%
      # 
      # add_trace(x = ~xnew,
      #           y = pnew[,"upr"],
      #           mode = "lines",
      #           line = list(color = input$EggerRibCol, width = 1.5),
    #           showlegend = F, name = "95% upper CI") %>%
    
    add_trace(x = ~Precision,
              y = ~Standardised.effect,
              mode = "markers",
              marker = list(color = input$EggerPtCol),
              showlegend = F, name = "Observed") %>%
      
      layout(
        xaxis = list(title = input$EggerXlab,
                     showgrid = F, showline = T,
                     linewidth = 1, rangemode = "tozero", zeroline = T,
                     zerolinecolor = toRGB("gray"), zerolinewidth = 1),
        yaxis = list(title = input$EggerYlab, showgrid = F, showline = T,
                     linewidth = 1, rangemode = "tozero", zeroline = T,
                     zerolinecolor = toRGB("gray"), zerolinewidth = 1),
        autosize = F,
        width = input$EggerWidth,
        height = input$EggerHeight
      )
  })
  
  })
