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
      ),
    plotOutput("plot_1"),
    verbatimTextOutput("debug")
    
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
output$debug <- renderPrint({
  cat("Selected:", input$knockout_mouse, "\n")
  d <- df_phynotype_statistical_scores_plot_1()
  cat("Rows:", nrow(d), "\n")
  str(d)
})

df_phynotype_statistical_scores <- data.frame(
    phenotype = rep(c("Phenotype1", "Phenotype2", "Phenotype3", "Phenotype4", "Phenotype5",
                      "Phenotype6", "Phenotype7", "Phenotype8", "Phenotype9", "Phenotype10"), times = 4),
    p_value = c(0.01,   0.45,  0.001, 0.72,  0.03,   0.58,  0.88,  0.21,  0.004, 0.67,
                0.60,   0.002, 0.38,  0.04,  0.91,   0.27,  0.009, 0.74,  0.52,  0.15,
                0.0001, 0.27,  0.55,  0.008, 0.66,   0.83,  0.36,  0.02,  0.47,  0.0006,
                0.33,   0.81,  0.02,  0.49,  0.0005, 0.007, 0.62,  0.95,  0.18,  0.41),
    knockout = rep(c("Mouse1", "Mouse2", "Mouse3", "Mouse4"), each = 10)
  )
  updateSelectInput(session, "knockout_mouse", choices = unique( df_phynotype_statistical_scores$knockout ) )
  
  df_phynotype_statistical_scores_plot_1 <- reactive (df_phynotype_statistical_scores %>% filter(. , knockout == input$knockout_mouse) )
  # df_phynotype_statistical_scores <- Here we should query the database to get a dataframe of p value + phenotype + knockout
  plot_1_title <- reactive(paste0( "The statistical scores of all recorded ", input$knockout_mouse  , " mouse phenotypes") )
    output$plot_1 <- renderPlot({
      # p-value - numerical continious. visualising over the phenotypes and gene knockout. Violin plot
      ggplot() + geom_col(data = df_phynotype_statistical_scores_plot_1(), aes(
        x = reorder(
          phenotype,
          p_value
        ) ,
        y = p_value
      ),
      colour = "black",
      fill ="coral"
      ) +
      theme_minimal() +
      labs(
        title = plot_1_title(),
        subtitle = "Smaller statistical value means higher correlation",
        x = "Mouse phenotype", # FILL IT
        y = " Statistical Score " # FILL IT
      )
    })
}

# Run the application 
runApp(shinyApp(ui = ui, server = server), launch.browser = TRUE)
