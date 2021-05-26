#'Explore your Data
#'
#'This function generates a shiny App that allows you to interactively explore your data.
#'
#'@param res A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes. Can also be a named list of results objects.
#'@param counts Normalized count matrix with rows as genes and columns as samples. Can also be a named list of count objects - ensure same order and names as the results object list.
#'@param cond Character vector indicating conditions belonging to each sample (same order as colnames(counts)). Can also be a named list of condition vectors - ensure same order and names as results and counts object lists.
#'@return A shiny app will be generated allowing you to explore your data interactively using BinfTools. Options to be selected:
#'        Dataset: the name of the dataset to look at (from your named list of objects). If no list is provided, defaults to 'data'.
#'        Control condition: The control condition for the current dataset (selected from cond)
#'        Absolute log2FoldChange Threshold: The log2 fold-change threshold for differentially expressed genes (defaults to log2(1.5))
#'        P-value column: Use adjusted p-value (padj) or raw p-value (pvalue)
#'        P-value threshold: The significance threshold for differentially expressed genes (default=0.05)
#'        Gene to colour/label: Type in gene name to colour/label on volcano plot, MA plot, and heatmap (must match gene symbols in res/counts objects)
#'        Scale to point sizes: If yes, points in volcano plot and MA plot will be scaled to gene expression levels or significance, respectively
#'        species: Indicate the species of the data for pathway enrichment analyses
#'        Pathway Enrichment Sources: Select data source for pathway enrichment analysis
#'        Print top 10 terms: sig= top 10 significant terms from pathway enrichment. enr= top 10 enriched terms from pathway enrichment
#'        Pathway Keyword (GO only): Type in a key word for targeted analysis. Must have no spaces. Term names containing this string (even as part of a word) will be included for targeted analysis
#'        Count plot method: ind= use values from individual samples, mean= use average from samples in each condition, geoMean= use geometric mean from samples in each condition, median= use median from samples in each condition
#'
#' @export

exploreData<-function(res, counts, cond){
  require(gpGeneSets, quietly = T)
  if(!is.list(res)){
    res<-list(data=res)
  }
  if(!is.list(counts)){
    counts<-list(data=counts)
  }
  if(!is.list(cond)){
    cond<-list(data=cond)
  }
  # Define UI for app that draws a histogram ----
  ui <- shiny::fluidPage(

    # App title ----
    shiny::titlePanel("Explore your Data"),

    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(

      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        shiny::selectInput(inputId="data",
                      label="Dataset",
                      choices=names(res),
                      selected=names(res)[1]),
        shiny::selectInput(inputId="control",
                           label="Select control condition:",
                           choices=unique(unlist(cond)),
                           selected=unique(unlist(cond))[1]),
        # Input: Numeric input for the log2FoldChange Threshold ----
        shiny::numericInput(inputId = "FC",
                     label = "Absolute log2FoldChange Threshold:",
                     min = 0,
                     max = NA,
                     value = log2(1.5)),
        # Radio buttons to determine if adjusted or raw p-value should be used
        shiny::radioButtons(inputId= "sig",
                     label="P-value column:",
                     choices=c("padj","pvalue"),
                     selected="padj"),

        #Input: Numeric input for the
        shiny::numericInput(inputId="p",
                     label = "P-value threshold:",
                     min=0,
                     max=1,
                     value=0.05),
        #Input: Text input for genes to colour
        shiny::textInput(inputId="colGene",
                  label="Gene to colour/label:",
                  value="",
                  placeholder = rownames(res[[1]])[1]),

        shiny::radioButtons(inputId="scale",
                     label="Scale point sizes?",
                     choices=c("Yes","No"),
                     selected="No"),
        shiny::radioButtons(inputId="species",
                            label="Species:",
                            choices=c("hsapiens","mmusculus", "dmelanogaster"),
                            selected="hsapiens"),
        shiny::selectInput(inputId="source",
                           label="Pathway Enrichment Sources:",
                           choices=c("GO","GO:BP","GO:CC","GO:MF",
                                     "KEGG","WP","REAC","TF","MIRNA"),
                           selected="GO"),
        shiny::radioButtons(inputId="GOprint",
                            label="Print top 10 terms:",
                            choices=c("sig","enr"),
                            selected="sig"),
        shiny::textInput(inputId="target",
                         label="Pathway Keyword (GO only):",
                         value="",
                         placeholder="Interleukin"),
        shiny::selectInput(inputId="method",
                           label="Count plot method:",
                           choices=c("ind","mean","geoMean","median"),
                           selected="ind")
      ),

      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::tabsetPanel(type="tabs",
                           shiny::tabPanel("DE Analysis",
                                           shiny::fluidRow(
                                             shiny::column(12, shiny::plotOutput(outputId="volcanoPlot")),
                                             shiny::column(12, shiny::plotOutput(outputId="MAPlot"))
                                           )),
                           shiny::tabPanel("Pathway Enrichment",
                                           shiny::fluidRow(
                                             shiny::column(12, shiny::plotOutput(outputId="GO_up")),
                                             shiny::column(12, shiny::plotOutput(outputId="GO_down"))
                                           )),
                           shiny::tabPanel("Targeted Pathay Analysis",
                                           shiny::fluidRow(
                                             shiny::column(12, shiny::plotOutput(outputId="gsva")),
                                             shiny::column(12, shiny::plotOutput(outputId="heatmap")),
                                             shiny::column(12, shiny::plotOutput(outputId="countPlot"))
                                           ))
        )
        # Output: Histogram ----
        # plotOutput(outputId = "volcanoPlot")
        # plotOutput(outputId= "MAPlot")
        #
      )
    )
  )


  server <- function(input, output) {

    get.DEG<-shiny::reactive({
      x<-res[[which(names(res) %in% input$data)]]
      #print(head(x))
      if(input$sig=="padj"){
        return(list(Up=rownames(subset(x, log2FoldChange >= input$FC & padj < input$p)),
             Down=rownames(subset(x, log2FoldChange <= -input$FC & padj < input$p))))
      }
      if(input$sig=="pvalue"){
        return(list(Up=rownames(subset(x, log2FoldChange >= input$FC & pvalue < input$p)),
             Down=rownames(subset(x, log2FoldChange <= -input$FC & pvalue < input$p))))
      }
    })
    GO_res<-shiny::reactive({
      x<-res[[which(names(res) %in% input$data)]]
      DEG<-get.DEG()
      #print(lapply(DEG, head))
      BinfTools::GO_GEM(DEG, species=input$species, bg=rownames(x),
                                source=input$source, prefix=paste(input$data,input$source, sep="_"),
                                pdf=F, fig=F, writeRes=F, writeGem=F, returnRes=T)
    })
    get.gmt<-shiny::reactive({
      x<-res[[which(names(res) %in% input$data)]]
      DEG<-get.DEG()
      #print(lapply(DEG, head))
      y<-BinfTools::GO_GEM(DEG, species=input$species, bg=rownames(x),
                        source=input$source, prefix=paste(input$data,input$source, sep="_"),
                        pdf=F, fig=F, writeRes=F, writeGem=F, returnRes=F, returnGost=T)
      if(input$species=="hsapiens"){
        return(c(BinfTools::customGMT(y$Up, input$target, gp_hs),
                 BinfTools::customGMT(y$Down, input$target, gp_hs)))
      }
      if(input$species=="mmusculus"){
        return(c(BinfTools::customGMT(y$Up, input$target, gp_mm),
                 BinfTools::customGMT(y$Down, input$target, gp_mm)))
      }
      if(input$species=="dmelanogaster"){
        return(c(BinfTools::customGMT(y$Up, input$target, gp_dm),
                 BinfTools::customGMT(y$Down, input$target, gp_dm)))
      }
    })
    output$volcanoPlot <- shiny::renderPlot({

      x    <- res
      if(is.list(res)){
        x<-res[[which(names(res) %in% input$data)]]
      }
      if(input$scale=="Yes"){
        scale<-T
      } else{
        scale<-F
      }

      if(input$sig=="padj"){
        BinfTools::volcanoPlot(x, title="Volcano", p=input$p, FC=input$FC,
                               lab=input$colGene, col=input$colGene, expScale=scale, returnDEG=F)
      }
      if(input$sig=="pvalue"){
        BinfTools::volcanoPlot(x, title="Volcano", pval=input$p, FC=input$FC,
                               lab=input$colGene, col=input$colGene, expScale=scale, returnDEG=F)
      }


    })
    output$MAPlot<-shiny::renderPlot({
      x<-res
      if(is.list(res)){
        x<-res[[which(names(res) %in% input$data)]]
      }
      if(input$scale=="Yes"){
        scale<-T
      } else{
        scale<-F
      }

      if(input$sig=="padj"){
        BinfTools::MA_Plot(x, "MA Plot", p=input$p, FC=input$FC, col=input$colGene,
                           lab=input$colGene, sigScale=scale)
      }
      if(input$sig=="pvalue"){
        BinfTools::MA_Plot(x, "MA Plot", pval=input$p, FC=input$FC, col=input$colGene,
                           lab=input$colGene, sigScale=scale)
      }

    })
    output$GO_up<-shiny::renderPlot({
      x<-GO_res()
      BinfTools:::GO_plot(x$Up, prefix=paste("Up",input$source, sep="_"),
                          ts=c(10,500), pdf=F, fig=T, print=input$GOprint)
    })
    output$GO_down<-shiny::renderPlot({
      x<-GO_res()
      BinfTools:::GO_plot(x$Down, prefix=paste("Down",input$source, sep="_"),
                          ts=c(10,500), pdf=F, fig=T, print=input$GOprint)
    })
    output$gsva<-shiny::renderPlot({
      x<-counts[[which(names(counts) %in% input$data)]]
      co<-cond[[which(names(cond) %in% input$data)]]
      gmt<-get.gmt()
      BinfTools::gsva_plot(t(scale(t(x))), gmt, method="ssgsea", condition=co, title=paste0("Pathway enrichment - ", input$target))
    })
    output$heatmap<-shiny::renderPlot({
      x<-counts[[which(names(counts) %in% input$data)]]
      co<-cond[[which(names(cond) %in% input$data)]]
      gmt<-unique(unlist(get.gmt()))
      BinfTools::zheat(genes=gmt, counts=x, conditions=co, con=input$control, title=paste0("Gene expression - ", input$target), labgenes=input$colGene)
    })
    output$countPlot<-shiny::renderPlot({
      x<-counts[[which(names(counts) %in% input$data)]]
      co<-cond[[which(names(cond) %in% input$data)]]
      gmt<-unique(unlist(get.gmt()))
      BinfTools::count_plot(x, genes=gmt, condition=co, title=paste0("Gene expression - ",input$target),
                            method=input$method)
    })
  }

  # Create Shiny app ----
  shiny::shinyApp(ui = ui, server = server)
}
