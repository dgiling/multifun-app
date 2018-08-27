
library(shiny)

# Functions
plot.simple <- function(data, x, y, xlab, ylab) {
  ggplot(data, aes_string(x=x, y=y)) +
    geom_point(size=3, alpha=0.2) +
    #labs(title=main) +
    scale_y_continuous(name=ylab) +
    scale_x_continuous(name=xlab) +
    theme_bw(base_size=15) +
    #stat_smooth(method="lm", colour="black", size=1, se=F) +
    theme(panel.grid = element_blank())
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# Define server logic
shinyServer(function(input, output, session) {
  shinyjs::runjs('toggleCodePosition();')

  observeEvent(input$gen.data, {
    updateTabsetPanel(session, "main",
                      selected = "Data simulation results")
  })
  
  ## Generate data table
  dt <- reactive ({   

    input$gen.data
    
    # Define GCD
    isolate({
      def <- defData(varname = "GCD", dist = "uniformInt", formula = paste0("1;",input$nlev))
    })
    
    # Step 1
    if(input$dir == "left") { 
      isolate({
      e1 <- renderText({ input$eqn1.l }, quoted=T)
      e1var <- input$eqn1var.l
      e2 <- renderText({ input$eqn2.l }, quoted=T)
      e2var <- input$eqn2var.l
      })
      def <- defData(def, dist="normal", varname = "M2", formula = e2(), variance = e2var) 
      def <- defData(def, dist="normal", varname = "M1", formula = e1(), variance = e1var)
    }
    if(input$dir == "right") { 
      isolate({
      e1 <- renderText({ input$eqn1.r }, quoted=T)
      e1var <- input$eqn1var.r
      e2 <- renderText({ input$eqn2.r }, quoted=T)
      e2var <- input$eqn2var.r
      })
      def <- defData(def, dist="normal", varname = "M1", formula = e1(), variance = e1var)
      def <- defData(def, dist="normal", varname = "M2", formula = e2(), variance = e2var) 
    }
  
    # Step 2
    isolate({
    e3 <- renderText({ input$eqn3 }, quoted=T)
    e3var <- input$eqn3var
    e4 <- renderText({ input$eqn4 }, quoted=T)
    e4var <- input$eqn4var
    e5 <- renderText({ input$eqn5 }, quoted=T)
    e5var <- input$eqn5var    
    e6 <- renderText({ input$eqn6 }, quoted=T)
    e6var <- input$eqn6var
    })
    
    def <- defData(def, dist="normal", varname = "F1", formula = e3(), variance = e3var)
    def <- defData(def, dist="normal", varname = "F2", formula = e4(), variance = e4var) 
    def <- defData(def, dist="normal", varname = "F3", formula = e5(), variance = e5var)
    def <- defData(def, dist="normal", varname = "F4", formula = e6(), variance = e6var) 
    
    # Generate data
    isolate({ sim.data <- genData(input$nobs, def) })

    # calculate EMF
    # 1. averaging
    vars <- c("F1", "F2", "F3", "F4")
    sim.data <- cbind(sim.data, getStdAndMeanFunctions(sim.data, vars)) # standardises between 0 and 1, takes average
    
    # 2. threshold
    vars.std <- c("F1.std", "F2.std", "F3.std", "F4.std")
    thresh <- getFuncsMaxed(sim.data, vars.std, threshmin=0.05, threshmax=0.99, prepend=c("id","GCD"))

    thresh25 <- thresh[thresh$thresholds=="0.25",c("id","GCD","funcMaxed")]
    names(thresh25)[3] <- "funcMaxed25"
    thresh50 <- thresh[thresh$thresholds==0.50,c("id","GCD","funcMaxed")]
    names(thresh50)[3] <- "funcMaxed50"
    thresh75 <- thresh[thresh$thresholds=="0.75",c("id","GCD","funcMaxed")]
    names(thresh75)[3] <- "funcMaxed75"
    
    sim.data <- merge(sim.data, thresh25, by=c("id","GCD"))
    sim.data <- merge(sim.data, thresh50, by=c("id","GCD"))
    sim.data <- merge(sim.data, thresh75, by=c("id","GCD"))
    
    sim.data

  })
    
    
  ## Generate plots
  # Output plot A: Step 1
  output$plotA <- renderPlot({
    
      plot1 <- plot.simple(data=dt(), x="GCD", y="M2", xlab="Global change driver", ylab="Mechanism 1 (M1)") #, main = "Equation 1")
      
      if(input$dir == "left") {  
       plot2 <- plot.simple(data=dt(), x="M2", y="M1", xlab="Mechanism 2 (M2)", ylab="Mechanism 1 (M1)") #, main = "Equation 1")
      }
      if(input$dir == "right") {  
       plot2 <- plot.simple(data=dt(), x="M1", y="M2", xlab="Mechanism 1 (M1)", ylab="Mechanism 2 (M2)") #, main = "Equation 2")
      }
       
      plot3 <- plot.simple(data=dt(), x="GCD", y="M2", xlab="Global change driver", ylab="Mechanism 2 (M2)") #, main = "Equation 2")
      plot_grid(plot1, plot2, plot3, ncol=3, align = 'h')

  })
  
  # Output plot B: Step 2, biotic effects
  output$plotB <- renderPlot({
  
      #vars.std <- c("F1.std", "F2.std", "F3.std", "F4.std")
      vars <- c("F1", "F2", "F3", "F4")
      
      f.long <- melt(dt(), id.vars=c("GCD","M2","M1"), measure.vars = vars)
      levels(f.long$variable) <- c("F1", "F2", "F3", "F4")
      
      ggplot(f.long, aes(x=M1, y=value)) +
        facet_wrap(~variable, ncol=4) +
        #geom_point(size=3, color = "#049C8E", alpha=0.3) + # biotic green
        geom_point(size=3, alpha=0.3) +
        scale_y_continuous(name="Function") +
        scale_x_continuous(name="Mechanism 1 (M1)") +
        theme_bw(base_size=15) +
        #stat_smooth(method="lm", colour="black", size=1, se=T) +
        theme(panel.grid = element_blank(), strip.background = element_blank())
      
  })
  
  # Output plot C: Step 2, abiotic effects
  output$plotC <- renderPlot({
    
    vars.std <- c("F1.std", "F2.std", "F3.std", "F4.std")
    vars <- c("F1", "F2", "F3", "F4")
    
    f.long <- melt(dt(), id.vars=c("GCD","M2","M1"), measure.vars = vars)
    levels(f.long$variable) <- c("F1", "F2", "F3", "F4")
    
    ggplot(f.long, aes(x=M2, y=value)) +
      facet_wrap(~variable, ncol=4) +
      #geom_point(size=3, color = "#FB9300", alpha=0.3) + # abiotic orange
      geom_point(size=3, alpha=0.3) +
      scale_y_continuous(name="Function") +
      scale_x_continuous(name="Mechanism 2 (M2)") +
      theme_bw(base_size=15) +
      #stat_smooth(method="lm", colour="black", size=1, se=F) +
      theme(panel.grid = element_blank(), strip.background = element_blank())
    
  })
  
  # plot D: correlation matrix
  output$plotD <- renderPlot({
    
    m <- cor(dt()[,c("F1", "F2", "F3", "F4")])
    p.mat <- cor.mtest(dt()[,c("F1", "F2", "F3", "F4")])
    
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot(m, method="color", col=col(200),  order="alphabet", number.cex=1.25,
             type="upper", outline=TRUE,
             addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=90, tl.cex = 1.5, #Text label color and rotation
             # Combine with significance
             p.mat = p.mat, sig.level = 0.05, insig = "blank", 
             # hide correlation coefficient on the principal diagonal
             diag=FALSE 
    )
    

  })
  
  # plotE: EMF
  output$plotE <- renderPlot({
    
    ggplot(dt(), aes(x=GCD, y=meanFunction)) +
      geom_point(size=3, alpha=0.2) +
      theme_bw(base_size=15) +
      scale_y_continuous(name="Mean function", limits=c(0,1)) +
      scale_x_continuous(name="Global change driver") +
      stat_smooth(method="lm", colour="black", size=1, se=T) +
      theme(panel.grid = element_blank())
    
  })
    
  
  output$plotF <- renderPlot({
    
    emf.long <- melt(dt(), id.vars="GCD", measure.vars = c("funcMaxed25", "funcMaxed50", "funcMaxed75"))
    levels(emf.long$variable) <- c("25% threshold", "50% threshold", "75% threshold")
    
    ggplot(emf.long, aes(x=GCD, y=value)) +
      facet_wrap(~variable) +
      geom_point(size=3, alpha=0.2) +
      scale_y_continuous(name="Number of functions") +
      scale_x_continuous(name="Global change driver") +
      theme_bw(base_size=15) +
      geom_smooth(method="glm" ,colour="black", size=1, se=T) +  # not an appropriate model
      theme(panel.grid = element_blank(), strip.background = element_blank())
    
  })
  
  # plotG: SEM
  output$plotG <- renderPlot({
    
        if(input$dir == "left") { 
        model.l <-  paste0("
                  # regressions
                  M1 ~ GCD + M2
                  M2 ~ GCD 
                  ",
                  input$SEMresponse, " ~ M1 + M2
                  ")
        fit <- sem(model.l, data=dt())
        }
        if(input$dir == "right") { 
          model.r <-  paste0("
                           # regressions
                           M1 ~ GCD
                           M2 ~ GCD  + M1
                           ",
                           input$SEMresponse, " ~ M2 + M1
                           ")
          fit <- sem(model.r, data=dt())
        }
        #summary(fit, standardized=TRUE)
        
        # plot SEM
        semPaths(fit, "std", residuals=F, edge.label.cex = 3, sizeMan = 20, sizeMan2=10,
                 #nodeLabels=c("M1","M2", "EMF", "GCD"), color=c("#049C8E","#FB9300","grey","grey"), # green/orange
                 nodeLabels=c("M1","M2", "EMF", "GCD"), color=c("grey","grey","grey","grey"),
                 layout=cbind(c(-1,1,0,0),c(0,0,-1,1)))
        
  })
})










