
# Packages
library(shiny)
library(shinyjs)
library(simstudy)
library(corrplot)
library(reshape2)
library(devtools)
library(multifunc)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(lavaan)
library(semPlot)
library(dplyr)
library(irlba)

# Define UI
shinyUI(
  fluidPage(
  useShinyjs(),
  # Title
  titlePanel("Generate hypotheses for global change effects on ecosystem multifunctionality"),
  
  fluidRow(
    column(3, style = "background-color:#F2F2F2;",

      column(12, 

        # input for Step 1
        br(),
        h5("Step 1"),
        radioButtons(inputId = "dir", "Direction of arrow a", choices =c("left", "right"), inline=T)
      ),
      
      column(8, 
        
        # input for Step 1
        # left
        conditionalPanel( condition = "input.dir == 'left'",
               textInput(inputId="eqn1.l", label="Eqn 1: mean M1", value="1*GCD + 0*M2"),
               textInput(inputId="eqn2.l", label="Eqn 2: mean M2", value="1*GCD")
        ),
      
        # right
        conditionalPanel( condition = "input.dir == 'right'",
                textInput(inputId="eqn1.r", label="Eqn 1: mean M1", value="1*GCD"),
                textInput(inputId="eqn2.r", label="Eqn 2: mean M2", value="1*GCD + 0*M1")
        ),
        
        # input for Step 2
        br(),
        h5("Step 2"),
        textInput(inputId="eqn3", label="Eqn 3: mean F1", value="1*M1 + 1*M2"),
        textInput(inputId="eqn4", label="Eqn 4: mean F2", value="1*M1 + 1*M2"),
        textInput(inputId="eqn5", label="Eqn 5: mean F3", value="1*M1 + 1*M2"),
        textInput(inputId="eqn6", label="Eqn 6: mean F4", value="1*M1 + 1*M2")
       
      ),
      
      column(4,  
               
           # Step 1 variances
           # left
           conditionalPanel( condition = "input.dir == 'left'",
                             textInput(inputId = "eqn1var.l", label="variance", value=8, width="75%"),
                             textInput(inputId = "eqn2var.l", label="variance", value=8, width="75%")
           ),
           
           # right
           conditionalPanel( condition = "input.dir == 'right'",
                             textInput(inputId = "eqn1var.r", label="variance", value=8, width="75%"),
                             textInput(inputId = "eqn2var.r", label="variance", value=8, width="75%")
           ),
        
           # Step 2 variances
           br(),
           br(),
           br(),
           textInput(inputId = "eqn3var", label="variance", value=8, width="75%"),
           textInput(inputId = "eqn4var", label="variance", value=8, width="75%"),
           textInput(inputId = "eqn5var", label="variance", value=8, width="75%"),
           textInput(inputId = "eqn6var", label="variance", value=8, width="75%")
           
      ),
      column(4,
        br(),
        textInput(inputId = "nobs", label="Number of observations", value = 50, width = "100%")
      ),
      column(4, 
        br(),
        numericInput(inputId = "nlev", label="Number of GCD levels", value = 5)
      ),
      column(4,
        br(),
        br(),
        br(),
        actionButton("gen.data", "Generate data")     
      ),
      
      column(12, 
             
        # input for Step 3
        br(),
        h5("Step 3"),
        radioButtons(inputId = "SEMresponse", "Calculation method for SEM", inline=F,
                     choiceNames =c("average", "25% threshold", "50% threshold", "75% threshold"),
                     choiceValues =c("meanFunction", "funcMaxed25", "funcMaxed50", "funcMaxed75") ),
        #actionButton("updateSEM", "Update SEM"),
        hr()
        
      )
    ),
    
    # main panel
    column(9,
           
      tabsetPanel(id="main",
          tabPanel("About",
        
             column(8,
              
               # update this in word doc  
               br(),  
               
               h4("Introduction"),
               helpText("This application is designed to assist researchers investigating the effects of global change on ecosystem multifunctionality (i.e. GC-EMF research) to develop mechanistic hypotheses. The application accompanies", a("Giling et al. (2018)", href="https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.14528"),". Please see the paper for a comprehensive description of the framework presented here."),
               helpText("The user provides hypothesized effect sizes for 1) the effect of a global change driver on biotic and/or abiotic mechanisms that drive functioning (e.g. species diversity or water availability), and 2) the influence of these mechanisms on multiple individual functions. The application then generates a simple dataset based on these hypotheses to explore the multifunctionality responses using both the averaging and threshold approaches. These simulations can be used to quickly visualize and disentangle how the interacting ecological and mathematical effects underlying the conceptual model influence the emergent GC-EMF relationship."),
               br(),
               
               h4("Instructions"),
               helpText("Hypotheses are entered in the sidebar on the left. Each equation represents one or two arrows from the conceptual model below. A zero (e.g. '0*M1' or just '0') indicates no effect. Interactive effects can also be added (e.g. '1*M1*M2'). The arrow 'a' can be set as pointing left or right (but not bidirectional). To exclude arrow a, simply do not include M1 or M2 as a predictor variable in Eqn. 1 or 2."),
               helpText("Click 'Generate data' after making hypotheses to run a simulation"),
               br()
                 
             ),
             
             column(12,
                           
               img(src='workflow_box_fig1b.png', align = "left")       
                     
             ),
             
             column(8, 
                    
               helpText("Figure 1. Conceptual framework demonstrating the suggested steps for building mechanistic hypotheses of global change effects on ecosystem multifunctionality. Step 1 entails defining the effect of a global change driver on mechanisms that may be important for ecosystem function (blue arrows). Two mechanisms are displayed here (M1 and M2), and they may influence each other (horizontal blue arrow a). At Step 2, the hypothesized additive or interactive effects of these mechanisms on individual functions (F1-F4) are defined (red arrows). Finally, Step 3 involves considering how the method of multifunctionality calculation affects the interpretation of emergent global change effects (green brace)."),     
                    
               br(),
               h4("Methods"),
               helpText("Data are generated according to the relationships specified by the user using the simstudy package (Goldfeld 2018). The correlation matrix displays positive (blue) and negative (red) correlations that are significant at the level p = 0.05. Multifunctionality metrics are calculated from the simulated data using the multifunc package (Byrnes et al. 2014; Byrnes 2017). Lines of model fit are modelled by stat_smooth in ggplot2 (Wickham 2009) with 'lm' and 'glm' for the averaging and threshold approach, respectively. Predictions are shown for visual aid only, and no checks are made on model fit or normality and dispersion assumptions. A structural equation model is fitted using lavaan (Rosseel 2012) and plotted with semPlot (Epskamp 2017) (again there are no checks on fit). All predictions will be invalidated when relationships are non-linear."),
               br(),
               
               h4("References"),
               helpText("Byrnes, J (2017). multifunc: Analysis of Ecological Drivers on Ecosystem Multifunctionality. R package version 0.8.0. https://github.com/jebyrnes/multifunc"),
               helpText("Byrnes, J. E. K. et al. Investigating the relationship between biodiversity and ecosystem multifunctionality: challenges and solutions. Methods in Ecology and Evolution 5, 111-124, doi:10.1111/2041-210X.12143 (2014)."),
               helpText("Epskamp, S and with contributions from Simon Stuber (2017). semPlot: Path Diagrams and Visual Analysis of Various SEM Packages' Output. R package version 1.1. https://CRAN.R-project.org/package=semPlot"),
               helpText("Giling D.P., Beaumelle L., Phillips H.R.P., Cesarz S., Eisenhauer N., Ferlian O., Gottschall F., Guerra C., Hines J., Sendek A., Siebert J., Thakur M.P. & Barnes A.D. (accepted) A niche for ecosystem multifunctionality in global change research. Global Change Biology."),
               helpText("Goldfeld, K (2018). simstudy: Simulation of Study Data. R package version 0.1.9. https://CRAN.R-project.org/package=simstudy"),
               helpText("Rosseel, Y (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/."),
               helpText("Wickham, H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009."),
               
               br()
               
             )  
          ),
         
          tabPanel("Data simulation results",
            
             column(12, align="left",
               br(),
               tags$b("Step 1: Simulated effect of global change on mechanisms (equations 1 and 2)"),
               br(),
               br(),
               plotOutput("plotA", width = "45%", height = 175), # Eqn 1/2
               hr()
             ),
             

             column(8, align="left",
               tags$b("Step 2: Simulated effect of mechanisms on individual functions (equations 3-6)", align="left"),
               plotOutput("plotB", width = "75%", height = 190), # Eqn 3-6
               plotOutput("plotC", width = "75%", height = 190),
               hr()
             ),
             
             column(4, align="center",
               h5("Correlation matrix for functions 1-4"),
               br(),
               br(),
               plotOutput("plotD", width = "100%", height = 300), # correlation matrix
               br(),               
               br(),
               hr()
             ), 
              
             column(9,
                    tags$b("Step 3: Effect of correlation structure and multifunctionality metric on GC-EMF relationship", align="center")
             ),
             column(3,
                    h5("Structrual Equation Model", align="center")
             ),
             
             column(3, align="center",
               h5("Averaging approach"),
               plotOutput("plotE", width = "90%", height = 250) # EMF
             ), 
             column(6, align="center",
               h5("Threshold approach"),
               plotOutput("plotF", width = "90%", height = 250) # EMF
             ), 
             column(3, align="center",
               plotOutput("plotG", width = "100%", height = 250) # SEM
             )
          )
        )
      )
    )
  )
)








