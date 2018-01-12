library(DT)
library(DNAcopy)
source("global.R")

# Setting the maximum file upload limit to 1 GB
options(shiny.maxRequestSize=1000*1024^2)

server = function(input, output, session){
  # First Bar
  roots <- getVolumes()
  shinyDirChoose(input, 'bam_dir', session=session, roots = roots)
  output$filepaths <- renderPrint({parseDirPath(roots, input$bam_dir)})
  
  df <- reactive({
    inFiles <- input$fileIn
    df <- data.frame()
    if (is.null(inFiles))
      return(NULL)
    for (i in seq_along(inFiles$datapath)) {
      tmp <- read.table(inFiles$datapath[i], header = T)  
      df <- rbind(df, tmp)
    }
    df
  })
  output$tbl <- DT::renderDataTable(
    df()
  )
  dataInput_bar1<- eventReactive(input$action_bar1,{
    inFiles <- input$fileIn
    df <- data.frame()
    if (is.null(inFiles))
      return(NULL)
    #Adding progress bar
    progress1 <- shiny::Progress$new()
    on.exit(progress1$close())
    progress1$set(message = "Data preprocessing", value = 0)
    for (i in seq_along(inFiles$datapath)) {
      progress1$inc(0.1, detail = paste("Reading bed file"))
      tmp <- read.table(inFiles$datapath[i], header = T)  
      df <- rbind(df, tmp)
    }
    df
    data_preproc(bed_file = df, chr = input$chr_selected, bam_dir = parseDirPath(roots, input$bam_dir), min_cov = input$min_cov, filename = input$cov_filename, progress = progress1)
  })
  ## Data Summary
  output$text_bam_dir<- renderText({
    paste(dataInput_bar1()[[1]])
  })
  output$text_bam_numbers<- renderText({
    paste(dataInput_bar1()[[2]])
  })
  output$text_bed_regions<- renderText({
    paste(dataInput_bar1()[[3]])
  })
  output$text_cov_summary<- renderText({
    paste(dataInput_bar1()[[4]])
  })
  output$text_cov_rdata<- renderText({
    paste(dataInput_bar1()[[5]])
  })
  
  # Second Bar
  #dataInput_bar2<- eventReactive(input$action_bar2,{
  #  source("functions/Player_Yr_Summary.R")
  #  Player_yr_summary(input$gender_bar2, input$yr_bar2, input$player_name_bar2)
  #})
  Cov_Mtx <- eventReactive(input$action_bar2,{
    #if(input$Cov_Mtx_bar2 == 1){
    #  load(paste(parseDirPath(roots, input$bam_dir), "/", input$cov_filename, ".RData", sep=""))
    #}
    #else{
      inFiles <- input$CovfileIn
      load(inFiles$datapath[1])
    #}
    #Adding progress bar
    progress2 <- shiny::Progress$new()
    on.exit(progress2$close())  
    progress2$set(message = "Running HCMM_CNVs", value = 0)
    HCMM_CNVs(Cov_matrix = matrix_adj_Cov_rm_duplicated, n_cluster = input$HC_n_clusters, bed_file_sorted = bed_file_sorted, sample_names = sample_names, filename = input$CBS_filename , progress = progress2)
  })
  
  output$test<- renderText({
    paste(Cov_Mtx()[[1]])
  })
  
  # Third Bar
  #observe({
  #eventReactive(input$CBSfileIn,{
  observeEvent(input$CBSfileIn,{
    #source("functions/Tourney_Search.R")
    inFiles <- input$CBSfileIn
    load(inFiles$datapath[1])
    label_samples<- names(CBS_all)
    updateSelectizeInput(session, "sample_bar3",
                         choices = label_samples, server = TRUE)
    #updateSelectInput(session, "tourney_bar3",
    #                  label = "Select Tourney", choices = label_tourneys)
  })
  
  #eventReactive(input$action_bar3,{
    #inFiles <- input$CovfileIn
    #load(inFiles$datapath[1])
    #output$plot1 <- renderPlot({
      #inFiles <- input$CBSfileIn
      #load(inFiles$datapath[1])
      #plot(CBS_all[[match(input$sample_bar3, names(CBS_all))]], xmaploc = T )
    #})
    #plot(CBS_all[[match(input$sample_bar3, names(CBS_all))]], xmaploc = T )
  #})
  #observeEvent(input$action_bar3,{
  #  output$plot1 <- renderPlot({
  #  input$action_bar3
  #  inFiles <- input$CBSfileIn
  #  load(inFiles$datapath[1])
  #  plot(CBS_all[[match(input$sample_bar3, names(CBS_all))]], xmaploc = T )
  #})
  #})
  
  data<- eventReactive(input$action_bar3, {
    inFiles <- input$CBSfileIn
    load(inFiles$datapath[1])
    plot_HCMMCNVs(result = CBS_all, sample_id = input$sample_bar3)
  })
  
  output$plot1 <- renderPlot({
    plot(data()[[1]], xmaploc = T )
  })
  
  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$download_fig <- downloadHandler(
    filename =  function() {
      paste("CNVs_", input$sample_bar3, ".png", sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      png(file)
      plot(data()[[1]], xmaploc = T ) # draw the plot
      dev.off()  # turn the device off
      
    } 
  )
  
}
