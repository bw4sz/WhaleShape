#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(xts)
library(leaflet)
library(dplyr)

date<-seq(as.Date("2015-01-01"), as.Date("2015-01-10"), by="day")
a<-xts(1:10,order.by=date)
df = data.frame(Lat = rnorm(1)+10, Long = rnorm(1),Id=a)

data_a<-data.frame(a)
data_a1<-data_a %>%  
  mutate("Lat" =as.numeric(df[1,1]),"Long"=as.numeric(df[2,1]),"Date"=rownames(data_a))
ui <- fluidPage(
  sliderInput("time", "date",min(date), 
              max(date),
              value = max(date),
              step=1,
              animate=T),
  leafletOutput("mymap")
)

server <- function(input, output, session) {
  points <- reactive({
    data_a1 %>% 
      filter(Date==input$time)
  })
  
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addMarkers(data = points(),popup=as.character(points()$a))
  })
}

shinyApp(ui, server)