library(shiny)
library(shinythemes)
library(shinyFiles)

shinythemes::themeSelector()
navbarPage(
  theme = shinytheme("cerulean"),
  "HCMMCNVs",
  # First bar: Title Search
  tabPanel("Data pre-processing",
           sidebarPanel(
             tags$div(tags$label(h4("1. Choose bam files directory"))),
             shinyDirButton("bam_dir", "Choose bam files", "Select directory of bam files"),
             tags$div(class="form-group shiny-input-container", 
                      tags$div(tags$label(h4("2. Bed file input"))),
                      tags$div(tags$label("Choose folder", class="btn btn-primary",
                                          tags$input(id = "fileIn", webkitfile = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                      #tags$label("No folder choosen", id = "noFile"),
                      tags$div(id="fileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                               tags$div(class="progress-bar")
                      )     
             ),
             selectInput("chr_selected", h4("3. Chromosome"),
                         choices = c(1:22), selectize = FALSE, selected = 19),
             numericInput("min_cov", label = h4("4. Minimum mean coverage"), value = 10),
             textInput("cov_filename", label = h4("5. Output file name"), value = "Test"),
             actionButton("action_bar1", "Run")
             #numericInput("number_clusters", label = h4("5. Number of clusters"), value = 3)
             #radioButtons("gender_bar1", label = h3("Gender"),
                          #choices = list("Men" = 1, "Women" = 2), 
                          #selected = 1),
             #sliderInput("yr_range_bar1", h3("Year Range:"),
                         #min = 1968, max = 2015, value = c(2009,2010)),
             #textInput("player_name_bar1", h3("Player's Name"), "Rafael Nadal"),
             #actionButton("action_bar1", "Update")
             #submitButton("Update")
           ),
           mainPanel(theme = "bootstrap.css",
             includeScript("./www/text.js"),
             tags$div(tags$label(h5("1. Bam files directory"))),
             verbatimTextOutput("text_bam_dir"),
             tags$div(tags$label(h5("2. Data Summary"))),
             textOutput("text_bam_numbers"),
             textOutput("text_bed_regions"),
             textOutput("text_cov_summary"),
             textOutput("text_cov_rdata"),
             tags$div(tags$label(h5("3. Selected bed file"))),
             tabPanel("Files table", dataTableOutput("tbl"))
             #textOutput("Message_bar1"),
             #tableOutput("Stat_Table_bar1")
             #h4("output$dir"),
             #verbatimTextOutput("dir"), br()
             #tableOutput("contents")
           )
  ),
  # Second bar: Hierarchical Clustering Mixture Model Copy Number Variants
  tabPanel("Run HCMMCNVs",
           sidebarPanel(
             #radioButtons("Cov_Mtx_bar2", label = h4("Coverage Matrix"),
            #              choices = list("Processed data" = 1, "Saved data" = 2), 
            #              selected = 1),
             tags$div(tags$label(h4("1. Load the coverage .RData"))),
             tags$div(class="form-group shiny-input-container", 
                      tags$div(tags$label("Choose Coverage RData", class="btn btn-primary",
                                          tags$input(id = "CovfileIn", webkitfile = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                      #tags$label("No folder choosen", id = "CovnoFile"),
                      tags$div(id="CovfileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                               tags$div(class="progress-bar")
                      )     
             ),
             numericInput("HC_n_clusters", label = h4("2. Hierarchical Clustering: number of clusters"), value = 3),
             tags$div(tags$label(h4("3. Cancer cell line only (optional)"))),
             radioButtons("radio_Ploidy", label = "Add ploidy input?",
                          choices = list("No" = 1, "Yes" = 2), 
                          selected = 1),
             fileInput("input_Ploidy_estimation", label = "Choose a file: "),
             textInput("CBS_filename", label = h4("4. Output file name"), value = "Test"),
             #selectInput('yr_bar2', h3('Year'),
             #           choices = c(1968:2015), selectize = FALSE),
             #textInput("player_name_bar2", h3("Player's Name"), "Rafael Nadal"),
             actionButton("action_bar2", "Run")
           ),
           mainPanel(
             textOutput("test")
             #textOutput("Message_bar2"),
             #tableOutput("Stat_Table_bar2")
           )
  ),
  # Third bar: Visualization
  tabPanel("Visulization",
           sidebarPanel(
             tags$div(tags$label(h4("1. Load the CBS result"))),
             tags$div(class="form-group shiny-input-container", 
                      tags$div(tags$label("Choose CBS results", class="btn btn-primary",
                                          tags$input(id = "CBSfileIn", webkitfile = TRUE, type = "file", style="display: none;", onchange="pressed()"))),
                      #tags$label("No folder choosen", id = "CBSnoFile"),
                      tags$div(id="CBSfileIn_progress", class="progress progress-striped active shiny-file-input-progress",
                               tags$div(class="progress-bar")
                      )     
             ),
             #selectInput("tourney_bar3", h3("Select Tourney"),
             #                choices = NULL),
             selectizeInput("sample_bar3", h4("Select sample"),
                            choices = NULL, multiple = F),
             actionButton("action_bar3", "Plot")
           ),
           mainPanel(
             plotOutput("plot1"),
             downloadButton(outputId = "download_fig", label = "Download the plot")
           )
  )
  
)
