
# packages ----------------------------------------------------------------

library(shiny)
library(shinythemes)
# library(shinyjs)
# library(Matrix)
# library(dplyr)
library(ggplot2)
# library(shinythemes)
# library(tidyr)
# library(BiocManager)
library(hdf5r)
# library(Seurat)
# library(HDF5Array)
# options(repos = BiocManager::repositories())
# BiocManager::install()
# library(DT)


# Load
# load(file = './data/appdata.RData') # This is faster by ~10%
source(file = 'helper.R')


# UI ----------------------------------------------------------------------

ui <- fluidPage(
  
  # jQuery chunk for ensuring DOM is fully loaded before user can interact
  tags$script(HTML(
    '$(document).ready(function () {
    $.getJSON("https://ipapi.co/json/", function (client) {
        Shiny.onInputChange("client", client);
    });
});'
  )),
  
  # Set layout colors/theme
  theme = shinytheme("yeti"),
  
  # App title
  titlePanel(
    fluidRow(
      column(
        width = 9, 
        div(
          h2(
            window_title, 
            style="display:inline;"
          ),
          tags$a(
            h3(
              title_link_text, 
              style="display:inline;"
            ),
            href = title_link_url, 
            target = "_blank"
          ),
          style='padding-left:10px;'
        )
      ),
      column(
        width = 3, 
        div(
          tags$a(
            h5(
              "Browser App", 
              style = "display:inline;color:black;vertical-align:middle;"
            ),
            href = "https://github.com/JaeLeeLab/sci_scRNAseq_portal", 
            target = "_blank"
          ),
          tags$a(
            img(
              height = 24, 
              width = 24, 
              src = "GitHub-Mark-64px.png", 
              style="display:inline;vertical-align:middle;"
            ),
            href = "https://github.com/JaeLeeLab/sci_scRNAseq_portal", 
            target="_blank"
          )
        ),
        align="right"
      )
    ),
    
    # This is to display on the headers of web browser.
    windowTitle = window_title
  ),
  
  # Create panel of tabs for different queries
  tabsetPanel(
    
    # About study panel
    tabPanel(
      title = 'About',
      br(),
      fluidRow(
        column(
          width = 7,
          br(),
          HTML(
            text = "<p style='font-size:15px'>This website accompanies the publication:<br><br>Lindsay M. Milich, James S. Choi, Christine Ryan, Susana R. Cerqueira, Sofia Benavides, Stephanie L. Yahn, Pantelis Tsoulfas, Jae K. Lee; Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord. <i>J Exp Med</i> 2 August 2021; 218 (8): e20210040. <a href='https://doi.org/10.1084/jem.20210040'>DOI: https://doi.org/10.1084/jem.20210040</a></p>"
          ),
          HTML(
            text = '<br><p><b>Click the "Gene Expression" tab to start exploring the data presented in the paper.</b></p>'
          ),
          HTML(
            text = '<br><h4>Background:</h4>The wound healing process that occurs after spinal cord injury (SCI) is critical for maintaining tissue homeostasis and limiting tissue damage, but eventually results in a scar-like environment that is not conducive to regeneration and repair. Within the first 7 days post-injury (dpi), multiple dynamic processes involving complex cellular heterogeneity occur, thus providing an opportunity for therpeutic intervention. This dataset provides a resource to better understand the cellular heterogeneity and interactions that occur in the SCI environment.'
          ),
          HTML(
            text = '<h4>Methods:</h4><p>C57BL/6J mice were subject to mid-thoracic contusion SCI and processed to generate single-cell suspensions. Cells were collected from uninjured, 1 dpi, 3 dpi, and 7dpi spinal cords. In total, 66,178 cells were sequenced. '
          ),
          # HTML(
          #   text = '<p><li><b>Gene expression:</b> type in a gene of interest to plot its expression in the UMAP low dimensional space.</li><p>'
          # ),
          # HTML(
          #   text = '<p><li><b>Gene expression time course:</b> type in a gene of interest to plot its expression in the UMAP low dimensional space, with cells from different injury time-points plotted separately.</li></p>'
          # ),
          # HTML(
          #   text = "<p><li><b>Cluster marker genes:</b> select a cluster of interest to display the genes that show preferential expression in that cluster over others. Clicking on a row in the table will show the corresponding expression plot. On the left, you can select the type of test, which controls how many clusters (all or 75%) have lower expression of the gene, compared to the selected cluster.</li></p>"
          # ),
          style='padding-left:10px;'
        ),
        column(
          width = 5,
          div(
            tags$a(
              img(
                height = 350, 
                width = 475, 
                src = "homepage_sci_umap_time.png", 
                style="display:inline;vertical-align:middle;"
              )
            )
          ),
          div(
            tableOutput(
              outputId = 'count_table'
            )
          )
        )
      ),
      hr(),
      fluidRow(
        h3(
          "Data availability", 
          style="display:inline;"
        ),
        br(),
        br(),
        column(
          width = 12,
          HTML(
            text = "<p>Code used to analyze the single-cell RNA-seq data from <i>Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord</i> are available on <a href='https://github.com/JaeLeeLab/sci_scRNAseq'>Github.</a> Direct download links can be found in the repo.</p><p>Raw data are available from the SRA (Sequence Read Archive) database under study <a href='https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP295673'>SRP295673</a>. Gene-count matrices are available from the Gene Expression Omnibus under accession <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162610'>GSE162610</a>. Relevant sample-level and cell-level metadata are available under the GEO accession metadata.</p>"
          )
        ),
        hr(),
        style='padding-left:10px;'
      ),
      style='padding-left:20px; padding-right:20px'
    ),
    # Tab to query features 
    tabPanel(
      title = 'Gene Expression',
      shinyjs::useShinyjs(),
      sidebarLayout(
        fluid = TRUE,
        # `choices = NULL` so that choices are set in server-side
        sidebarPanel(
          tags$style(type = 'text/css', ".selectize-input { font-size: 16px; line-height: 22px;} .selectize-dropdown { font-size: 16px; line-height: 22px; }"),
          width = 2,
          selectizeInput(
            inputId = "selected_dataset",
            choices = NULL,
            label = 'Select dataset:'
          ),
          selectizeInput(
            inputId = "selected_feature", 
            choices = NULL, 
            options = list(placeholder = 'Select a gene'),
            label = 'Select gene:'
          ),
          selectizeInput(
            inputId = "selected_groupby",
            choices = NULL,
            label = 'Group cells by:'
          ),
          checkboxInput(
            inputId = 'draw_labels',
            label = 'Overlay group labels',
            value = TRUE
          ),
          shinyjs::hidden(
            textInput(
              inputId = "hidden_selected_feature", 
              label = "Hidden Selected Gene", 
              value = NULL
            )
          ),
          shinyjs::hidden(
            textAreaInput(
              inputId = "hidden_feature_list", 
              label = NULL, 
              value = NULL
            )
          ),
          shinyjs::hidden(
            textInput(
              inputId = "hidden_selected_cluster",
              label = "Hidden Selected Cluster",
              value = NULL
            )
          ),
          HTML(
            text = "<p style='font-size:14px'><b>How to use:</b><br>To get started, selected a dataset from the \"Select dataset\" drop-down menu. <ul style=\"padding-left:10px\"><li>All_SCI: all cells from study.</li><li>Myeloid: macrophages, microglia, etc. (as in Figure 2 of paper).</li><li>Vascular: endothelial cells, fibroblasts, etc. (as in Figure 6 of paper).</li><li>Macroglia: astrocytes, oligodendrocytes, etc. (as in Figure 8 of paper).</li></ul></p><p>Query genes from the \"Select Gene\" drop-down menu. Alternatively, start typing your gene of interest for matching items.</p><p>Cells can be grouped and colored by compartment, celltype, subtype, and other categorial meta-data, from the \"Group cells by\" menu.</p><p>To zoom in, click and drag to create a region of interest in either of the top two panels. Double-click the region to zoom. Double-click again to reset zoom.</p>"
          )
        ),
        mainPanel(
          fluid = TRUE, 
          width = 10,
          fluidRow(
            width = 12,
            column(
              width = 5,
              plotOutput(
                outputId = "dimplot",
                width = "100%", 
                height = "500px",
                click = 'dimplot_click',
                dblclick = 'dimplot_dblclick',
                brush = brushOpts(
                  id = 'umap_brush',
                  resetOnNew = TRUE
                )
              )
            ),
            column(
              width = 5,
              plotOutput(
                outputId = "featureplot",
                width = "auto", 
                height = "500px",
                click = 'splitfeatureplot_click',
                dblclick = 'featureplot_dblclick',
                brush = brushOpts(
                  id = 'umap_brush',
                  resetOnNew = TRUE
                )
              )
            ),
            column(
              width = 2,
              plotOutput(
                outputId = 'dimplotlegend',
                width = 'auto',
                height = '500px'
              )
            )
          ),
          fluidRow(
            width = 12,
            plotOutput(
              outputId = 'splitfeatureplot',
              width = 'auto',
              height = '350px',
              click = 'splitfeatureplot_click',
              dblclick = 'splitfeatureplot_dblclick',
              brush = brushOpts(
                id = 'umap_brush',
                resetOnNew = TRUE
              )
            )
          ),
          fluidRow(
            column(
              width = 3,
              # DTOutput(outputId = 'groupby_feature_table')
            ),
            column(
              width = 6,
              plotOutput(
                outputId = 'featuredotplot',
                hover = 'dotplot_hover',
                height = '550px'
              )
            ),
            column(
              width = 3,
              # DTOutput(outputId = 'groupby_feature_table')
            ),
            # column(
            #   width = 7,
            #   DTOutput(outputId = 'groupby_feature_table')
            # )
          )
        )
      )
    )
  )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  #Logging
  observeEvent(
    eventExpr = {
      input$client 
    },
    handlerExpr = {
      logging::loginfo("New client with ip: %s", input$client$ip)
    },
    ignoreNULL = TRUE, 
    ignoreInit = FALSE
  )
  
  # Cell count table
  output$count_table <- renderTable(
    expr = {
      cell_counts
    },
    bordered = TRUE,
    width = '100%'
  )
  
  # Initialize dataset
  updateSelectizeInput(
    session = session, 
    label = 'Select dataset:',
    inputId = "selected_dataset",
    choices = names(dataset_dict),
    selected = 'All_SCI'
  )
  
  # Initialize cell groupby
  updateSelectizeInput(
    session = session,
    inputId = "selected_groupby", 
    label = "Group cells by:",
    choices = cell_groupings,
    selected = 'Celltype'
  )
  
  # Initialize selected feature
  updateSelectizeInput(
    session = session,
    inputId = "selected_feature", 
    label = "Select Gene:",
    choices = all_features,
    server = TRUE,
    selected = 'Cx3cr1'
  )
  
  # Upon change to selected_dataset, take subset of expression matrix
  dataset_value <- eventReactive(
    eventExpr = {
      input$selected_dataset
    },
    valueExpr = {
      dataset_value <- dataset_dict[input$selected_dataset]
      return(dataset_value)
    },
    ignoreInit = FALSE,
    ignoreNULL = TRUE
  )
  
  
  # Upon change to selected_dataset, change point size according to number of 
  # points displayed.
  dataset_ptsize <- eventReactive(
    eventExpr = {
      input$selected_dataset
    },
    valueExpr = {
      return(pt_size_list[[dataset_value()]])
    }
  )
  
  dimplotlegend_labelsize <- eventReactive(
    eventExpr = {
      input$selected_groupby
    },
    valueExpr = {
      n <- length(unique(obs_sci[[input$selected_groupby]]))
      labelsize <- (n^2 + n) / (n^2 - n)
      return(labelsize)
    }
  )
  
  
  # Zooming on umaps
  ranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(
    eventExpr = {
      c(input$dimplot_dblclick, input$featureplot_dblclick, input$splitfeatureplot_dblclick)
    }, 
    handlerExpr = {
      gc(verbose = FALSE)
      brush <- input$umap_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    }
  )
  
  # Reset brush ROI on zoom
  observeEvent(
    eventExpr = {
      c(input$selected_dataset)
    }, 
    handlerExpr = {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  )
  
  # Dimplot to show clusters or other categorical metadata
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_groupby)
    },
    handlerExpr = {
      tmp_group <- req(input$selected_groupby)
      logging::loginfo("loaded dataset %s with dimplot cells group by %s.", 
                       dataset_value(), tmp_group)
      output$dimplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          DrawDimPlot(
            dataset = dataset_value(),
            groupby = tmp_group,
            ptsize = dataset_ptsize(),
            xranges = ranges$x,
            yranges = ranges$y,
            draw_labels = input$draw_labels
          )
        }
      )
    }
  )
  
  # Feature plot (to display gene expression)
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_feature)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature)
      output$featureplot <- renderPlot(
        expr = {
          logging::loginfo("loaded dataset %s with featureplot feature %s.", 
                           dataset_value(), tmp_feature)
          gc(verbose = FALSE)
          DrawFeaturePlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            ptsize = dataset_ptsize(),
            xranges = ranges$x,
            yranges = ranges$y
          )
        }
      )
    }
  )
  
  # Render dimplot legend, re-render upon change to `input$selected_groupby` or
  # `input$selected_dataset`.
  observeEvent(
    eventExpr = {
      c(input$selected_groupby, input$selected_dataset)
    },
    handlerExpr = {
      tmp_groupby <- req(input$selected_groupby)
      output$dimplotlegend <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          DrawDimPlotLegend(
            dataset = dataset_value(),
            groupby = tmp_groupby,
            ptsize = 5*dataset_ptsize(),
            labelsize = dimplotlegend_labelsize()
          )
        }
      )
    }
  )
  
  # Split feature plot - separate plot per injury time-point
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_feature)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature)
      output$splitfeatureplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          DrawSplitFeaturePlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            ptsize = dataset_ptsize(),
            xranges = ranges$x,
            yranges = ranges$y
          )
        }
      )
    }
  )
  
  # Feature dot plot with dot rows split by `input$selected_groupby`
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_feature, input$selected_groupby)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature)
      tmp_groupby <- req(input$selected_groupby)
      output$featuredotplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          DrawFeatureDotPlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            groupby = tmp_groupby
          )
        }
      )
    }
  )
}

shinyApp(ui = ui, server = server)

# To deploy, run the following two lines:
# Current working diretory should contain the project directory.
# setwd('..')
# rsconnect::deployApp(appDir = 'sci_scRNAseq_portal/', appName = 'sci_singlecell', account = 'jaeleelab')
# rsconnect::accounts()
# rsconnect::accountInfo()


# To allow bioconductor packages to be sourced. Run the following in your 
# current session before pushing.
# library(BiocManager)
# options(repos = BiocManager::repositories())
