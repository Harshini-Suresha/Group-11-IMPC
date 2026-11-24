install.packages("shiny", quietly = TRUE)
install.packages("DBI", quietly = TRUE)
install.packages("RMariaDB", quietly = TRUE)
install.packages("ggplot2", quietly = TRUE)
install.packages("tidyverse", quietly = TRUE)

library(shiny)
library(DBI)
library(RMariaDB)
library(ggplot2)
library(tidyverse)


# There are 3 visualisations that they require.
#(1) They want the ability to select a particular knockout mouse and visualise
# the statistical scores of all phenotypes tested. The plot should intuitively
# show the phenotypes that are significantly affected by the gene knockout.
# (2) The collaborator also wishes to visualise the
#statistical scores of all knockout mice for a selected phenotype. 
#(3) Lastly, they would like
#to visualise clusters of genes with similar phenotype scores. 



# Define UI for application
ui <- fluidPage(

    # Application title
    titlePanel("Placeholder"),

    # For (0) - they want to select a particular knockout mouse 
    selectInput( "knockout_mouse",
                  "Select particular knockout mouse:",
                  choices = NULL,
                  multiple = FALSE
                  )
    
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  df_phynotype_statistical_scores <- data.frame(
    phenotype = c("Phenotype1", "Phenotype2", "Phenotype3"),
    p_value = c(0.01, 0.05, 0.001),
    knockout = c("Mouse1", "Mouse2", "Mouse3")
  )
  updateSelectInput(session, "knockout_mouse", choices = unique( df_phynotype_statistical_scores$knockout ) )
  
  df_phynotype_statistical_scores_plot_0 <- reactive (df_phynotype_statistical_scores %>% dplyr::select(. , knockout = input$knockout_mouse) )
  # df_phynotype_statistical_scores <- Here we should query the database to get a dataframe of p value + phenotype + knockout
    output$plot_0 <- renderPlot({
      # p-value - numerical continious. visualising over the phenotypes and gene knockout. Violin plot
      ggplot() + geom_violin(data = df_phynotype_statistical_scores, aes(
        x = phenotype ,
        y = p_value
      ))
    })
}

# Run the application 
runApp(shinyApp(ui = ui, server = server), launch.browser = TRUE)
