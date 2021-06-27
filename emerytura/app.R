rm(list = ls())

library(shiny)
library(dummies)
library(MASS)
library(oglmx)
library(pscl)
library(ggplot2)
library(sandwich)
library(msm)
library(psych)
library('plot.matrix')
library(haven)


ui <- fluidPage(
  
  titlePanel("Jak dĹ‚ugo nalezy oszczedzac na godziwa emeryture?"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput('year','Rok urodzenia',min=1958,max=2019,value=2000),
      sliderInput('stopa','i - stopa procentowa',min=-0.05,max=0.05,value=-0.01,step=0.001),
      sliderInput('amin','minimalny wiek wejscia na rynek pracy',min=20,max=49,value=25,step=1),
      sliderInput('rmax','maksymalny wiek przejscia na emeryture',min=50,max=100,value=65,step=1),
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Opis Modeli",   withMathJax(),
                 h2('Plan z deterministycznym okresem akumulacji (Plan1)'),
                 br(),
                 h4('Modelowany wzorem: $$C \\cdot \\ddot{a}_{\\overline{r-a \\mid}} \\cdot (1+i)^{r-a} = E\\cdot \\ddot{a}_r  $$'),
                 br(),
                 h4('Gdzie:'),
                 br(),
                 h5('\\(a\\) to wiek wejĹ›cia na rynek pracy (wejĹ›cia do planu emerytalnego),'),
                 h5('\\(r\\) to wiek przejĹ›cia na emeryturÄ™,'),
                 h5('\\(E\\) to wypĹ‚acana emerytura (wyraĹĽona w procencie otrzymywanego wynagrodzenia),'),
                 h5('\\(C\\) to skĹ‚adka emerytalna (wyraĹĽona w procencie otrzymywanego wynagrodzenia),'),
                 h5('\\(\\ddot{a}_{\\overline{r-a \\mid}}\\) to renta finansowa dla fazy kumulacji kapitaĹ‚u,'),
                 h5('\\(\\ddot{a}_r\\) to renta Ë™ĹĽyciowa dla fazy emerytalnej.'),
                 hr(),
                 h2('Plan z ubezpieczeniowym okresem akumulacji (Plan2)'),
                 br(),
                 h4('Modelowany wzorem: $$\\frac{C \\cdot \\ddot{a}_{a:\\overline{r-a \\mid}} \\cdot (1+i)^{r-a}}{(r-a)p_a}= E\\cdot \\ddot{a}_r  $$'),
                 br(),
                 h4('Gdzie:'),
                 br(),
                 h5('\\(\\ddot{a}_{a:\\overline{r-a \\mid}}\\) to terminowa renta Ë™ĹĽyciowa')
),
        tabPanel("Plan 1",   withMathJax(),h2(textOutput("rok1")),
                 helpText('Wykres w formie map ciepła przedstawiajace jak zmieni sie stopa zastapienia
(skala nad wykresem 1 = 100% wynagrodzenia) w zaleznosci od parametrów \\(a\\), \\(r\\) i stopy procentowej przy
ustalonym \\(C\\)=15%.'),plotOutput("plotCall",width="100%")), 
        tabPanel("Plan 2",withMathJax(),h2(textOutput("rok2")),
                 helpText('Wykres w formie map ciepła przedstawiajace jak zmieni sie stopa zastapienia
(skala nad wykresem 1 = 100% wynagrodzenia) w zaleznosci od parametrów \\(a\\), \\(r\\) i stopy procentowej przy
ustalonym \\(C\\)=15%.'), plotOutput("plotPut",width="100%")) 
      )
    )
  )  
)



server <- function(input, output) {
  
  #Generate Black-Scholes values

  
  Plan1 = function(yr, s, min ,max) {
    MUP1 <- read_dta("MUP1.dta")
    MUP1 = subset(MUP1,MUP1$year  ==yr)
    MUP1$px = 1- MUP1$qx
    
    
    A1 <- matrix(0,nrow = 100, ncol = 100)
    B1 <- matrix(0,nrow = 100, ncol = 100)
    C1 <- matrix(0,nrow = 100, ncol = 100)
    
    r <- c(min:max)
    for (x in r)
    {
      a <- c(20:x)
      for (y in a) {
        z <- x-1
        i <- y:z 
        i <- -i+y
        A1[x,y]  <- sum((1+s)^i)
        j <- x:99 
        u=0
        for (t in j) {
          u <- u+prod(MUP1$px[(x+1):(t+1)])*((1+s)^(-t+x-1))
        }
        C1[x,y]  <- u+1
        B1[x,y]  <- 0.15*A[x,y]*(1+s)^(x-y)/C[x,y]
      }
    }
    
    D1=B1[min:max,min:max]
    colnames(D1)<-c(min:max)
    rownames(D1)<-c(min:max)
    res = D1
  }
  
  #Call plot
  output$plotCall <- renderPlot({
    yr=input$year
    s=input$stopa 
    min=input$amin
    max=input$rmax
    plot(Plan1(yr, s,min,max), border=NA, breaks=30, col=heat.colors(brk), 
         key=list(side=3,  font=2, cex.axis=0.75), fmt.key="%.2f", 
         polygon.key=NULL, axis.key=NULL, spacing.key=c(3,2,2),main="" , xlab="a - wiek wejscia na rynek pracy",ylab="r - wiek przejscia na emeryture")
  }, height = 350, width = 600)
  
  output$rok1<-renderText(paste("Plan z deterministycznym okresem akumulacji dla osoby urodzonej w ", input$year))
  output$rok2<-renderText(paste("Plan z ubezpieczeniowym okresem akumulacji dla osoby urodzonej w ", input$year))
  
  Plan2 = function(yr, s, min ,max) {
    MUP1 <- read_dta("MUP1.dta")
    MUP1 = subset(MUP1,MUP1$year  ==yr)
    MUP1$px = 1- MUP1$qx
    
    A <- matrix(0,nrow = 100, ncol = 100)
    B <- matrix(0,nrow = 100, ncol = 100)
    C <- matrix(0,nrow = 100, ncol = 100)
    V <- matrix(0,nrow = 100, ncol = 100)
    
    r <- c(min:max)
    
    for (x in r)
    {
      a <- c(20:x)
      for (y in a) {
        n=0
        k <- (y+1):(x+1) 
        for (g in k) {
          n <- n+prod(MUP1$px[(y+1):(g+1)])*((1+s)^(-g+y-1))
        }
        A[x,y]  <- n+1
        V[x,y]=prod(MUP1$px[(y+1):(x)])
        j <- x:99 
        u=0
        for (t in j) {
          u <- u+prod(MUP1$px[(x+1):(t+1)])*((1+s)^(-t+x-1))
        }
        C[x,y]  <- u+1
        B[x,y]  <- 0.15*A[x,y]*(1+s)^(x-y)/C[x,y]/V[x,y]
      }
    }
    
    D2=B[min:max,min:max]
    colnames(D2)<-c(min:max)
    rownames(D2)<-c(min:max)
    res = D2
  }
  

  #Call plot
  output$plotPut <- renderPlot({
    yr=input$year
    s=input$stopa 
    min=input$amin
    max=input$rmax
    plot(Plan2(yr, s,min,max), border=NA, breaks=30, col=heat.colors(brk), 
         key=list(side=3,  font=2, cex.axis=0.75), fmt.key="%.2f", 
         polygon.key=NULL, axis.key=NULL, spacing.key=c(3,2,2),main="", xlab="a - wiek wejscia na rynek pracy",ylab="r - wiek przejscia na emeryture")
  }, height = 350, width = 600)

  
  

  
}

shinyApp(ui, server)