################################################################################

### Shiny app - PPP UKB

################################################################################

# Developed by Danni A Gadd (https://github.com/DanniGadd)
# Email - danni.gadd@ed.ac.uk

# Load packages
library(shiny)
library(tidyverse)
library(readxl)
library(thematic)
library(ggplot2)
library(gridExtra)
library(networkD3)
library(igraph)
library(patchwork)
library(shinybusy)

# NOTE: Please update this path to point to the data folder
base_path <- "/srv/shiny-server/"

# Set annnotations
Text1 <- paste(" Circular points indicate associations with P < 3.1x10-6 (Bonferroni threshold), whereas triangles show no association. Failures in the Cox PH assumption (Schoenfeld P < 0.05) at the local protein level are shown in red.")

##############################################################################################

### SHINY APP CODE

### DEFINE UI

ui <- fluidPage(
    tags$head(tags$style(HTML("
        .selectize-input, .selectize-dropdown {
          font-size: 75%;
        }
        "))),

    theme = bslib::bs_theme(bootswatch = "minty"),

    headerPanel("Blood protein levels predict leading incident diseases and mortality in the UK Biobank"),
    add_busy_bar(color = "#33B9FF"),
    tabsetPanel(

        tabPanel("Cox associations",
            sidebarLayout(
                sidebarPanel(
                    width = 4,
                    uiOutput("protein_input"),
                    uiOutput("disease_input")
                ),

                mainPanel(
                    fluidRow(
                        plotOutput("plot2")
                    ),
                    fluidRow(
                        verbatimTextOutput("dfStr")
                    )
                )
            )),

        tabPanel("Network view",
                 sidebarLayout(
                    sidebarPanel(
                        width = 4,
                        uiOutput("protein_cols_input")
                     ),
                     mainPanel(
                         forceNetworkOutput(outputId = "plot1")
                     )
                 )
        )
    )
)

### DEFINE SERVER LOGIC

server <- function(input, output, session) {
    bind <- reactive({
        df <- read_csv(paste0(base_path, "latest_data/bind_narm.csv"))
        return(df)
    })

    assocs <- reactive({
        df <- read_csv(paste0(base_path, "latest_data/assocs.csv"))
        return(df)
    })

    assoc_proteins <- reactive({sort(unique(assocs()$Predictor))})
    output$protein_cols_input <- renderUI(
        selectInput(
            "Columns_protein",
            "Protein",
            selected = "GDF15.Q99988.OID20251.v1",
            choices = assoc_proteins(),
            multiple = TRUE
        )
    )

    protein_names <- reactive({sort(unique(bind()$Predictor))})
    output$protein_input <- renderUI(
        selectInput("p", "Protein", choices = protein_names())
    )

    disease_names <- reactive({sort(unique(bind()$Dis_name))})
    output$disease_input <- renderUI(
        selectInput("dis", "Disease", choices = disease_names())
    )

    ### TAB ONE - network view
    subset <- reactive({
        req(input$Columns_protein)
        table <- assocs()[which(assocs()$Predictor %in% input$Columns_protein),]
        table$Importance = ntile(table$HR, 6)
        table <- table[,which(colnames(table) %in% c('Outcome', 'Protein', 'Importance'))]
        return(table)
    })

    nodes <- reactive({
        nod1 <- subset()[,which(colnames(subset()) == 'Outcome')]
        nod1 <- as.data.frame(unique(nod1))
        names(nod1) <- 'name'
        nod1$group <- 'Outcome'
        nodp <- as.data.frame(unique(subset()$Protein))
        names(nodp) <- 'name'
        nodp$group <- 'Protein'
        nods <- rbind(nod1, nodp)
        nods$size <- 1
        return(nods)
    })

    links <- reactive({
        lnks <- subset()
        names(lnks) <- c('source', 'target', 'value')
        return(lnks)
    })

    net <- reactive({
        network <- graph_from_data_frame(d=links(), vertices=nodes(), directed=F)
        i_one_D3 <- igraph_to_networkD3(network, group = as.factor(V(network)$group), what = "both")
        return(i_one_D3)
    })

    output$plot1 <- renderForceNetwork(forceNetwork(
        Links = net()$links, Nodes = net()$nodes,
        Source = 'source', Value  = "value", Target = 'target', NodeID = 'name',
        Group = 'group', zoom = TRUE, opacity = 1, colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")))


    #####################################################################

    ### TAB TWO - year of follow up

    # Subset to chosen protein and disease
    filtered_data <- reactive({
        req(
            input$dis,
            input$p
        )
        test <- bind() %>%
            filter(
                Dis_name == input$dis,
                Predictor == input$p
            )
        return(test)
    })

    # Plot
    output$plot2 <- renderPlot({
        prot_name  <- as.character(unique(filtered_data()$Prot_name))
        p1 <- ggplot(filtered_data(),aes(Iteration, Hazard.Ratio))+
            geom_point(aes(shape = Association), fill = filtered_data()$col_variable, colour = filtered_data()$col_variable, size = 5,
                       shape = filtered_data()$shape_variable)+
            geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = filtered_data()$col_variable,
                          position = position_dodge(0.5), width = 0, linewidth = 0.75) +
            ylab("Hazard Ratio [95% CI]")+ xlab ("")+ theme_classic() +
            geom_hline(yintercept = 1, linetype = "dotted")+
            theme(title = element_text(colour = "aquamarine3"), axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "bottom",
                  axis.text.y = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
            theme_classic() + ggtitle(paste0(prot_name, ' - ', input$dis)) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())

        plot_data <- filtered_data() %>%
            dplyr::select('Iteration', 'Cases', 'Controls') %>%
            mutate(Year = Iteration,
                   Cases = Cases,
                   Controls = Controls) %>%
            pivot_longer(c(Year, Cases, Controls), names_to = "layer", values_to = "label")

        plot_data$layer <- factor(plot_data$layer, levels = c("Controls", "Cases", "Year"))


        p2 <-  ggplot(plot_data, aes(x = Iteration)) +
            geom_text(aes(y = factor(layer, c("Controls", "Cases", "Year")), label = label)) +
            labs(y = "", x = NULL) +
            theme_minimal() +
            theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
                  panel.grid = element_blank(), strip.text = element_blank())

        p1 / p2 +  plot_layout(heights = c(4, 1))


    }, res = 96, width = 900)

    # Print coption for annotations
    output$dfStr <- renderPrint({
        cat(Text1)
    })
}

# Run the application
shinyApp(ui = ui, server = server)




