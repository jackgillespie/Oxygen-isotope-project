## Version 0.1 Stable isotope data reduction ##
#install.packages("timetk")
.libPaths("C:/Users/Jack/Documents/R/lib")
install.packages("timetk")
#.libPaths("D:/R/lib2")
library(dplyr)
library(tidyverse)
library(stringr)
library(miniUI)
library(shiny)
library(timetk)
### test loop for entire folder ###
session_label <- "ExampleOdata"

#setwd(paste0("D:/R/O_isotope_UI/", session_label))
setwd(paste0("C:/Users/Jack/Documents/R/", session_label))

filename <- list.files(pattern = "*.asc")

combined_df <- c()

for (i in 1:length(filename)) {
  file_contents <- readLines(filename[i])
  
  # Parse acquisition filename
  line_number <- grep(file_contents,pattern="ACQUISITION FILE NAME")
  position <- str_locate(file_contents[line_number],"\\t")
  acquisition_filename <-  str_sub(file_contents[line_number],position[1]+1)
  spotname <- sub("^.+\\\\", "", acquisition_filename)
  
  # Parse session date
  line_number <- grep(file_contents,pattern="CAMECA")
  session_date <- as.character(str_match(file_contents[line_number],"[:digit:].*"))
  # Parse session time
  line_number <- grep(file_contents,pattern="CMCA")
  session_time <- as.character(str_match(file_contents[line_number],"[:digit:].*"))
  session_time <- format(strptime(session_time, "%I:%M %p"), format="%H:%M:%S")
  #combined for datetime
  session_datetime <- paste(session_date, session_time, sep = " ")

  session_datetime <- as.POSIXct(strptime(session_datetime, "%d/%m/%Y %H:%M:%S"))
  
  # Parse position coordinates
  line_number <- grep(file_contents,pattern="X POSITION")
  temp<-str_split(file_contents[line_number],"\\t",simplify=TRUE);
  x_position<-as.integer(temp[2]);
  y_position<-as.integer(temp[4]);
  
  
  
  # Parse the number of blocks
  line_number <- grep(file_contents,pattern="CUMULATED RESULTS")
  position <- str_locate(file_contents[line_number],"\\t[:digit:].*\\t")
  number_of_blocks <- as.integer(str_sub(file_contents[line_number],position[1]+1,position[2]-1))
  
  # Parse the raw data table
  line_number <- grep(file_contents,pattern="RAW DATA")
  column_names <- cbind(str_split(str_squish(file_contents[line_number+2])," ",simplify=TRUE),str_split(str_squish(file_contents[line_number+4])," ",simplify=TRUE))
  raw_data <- read.table(text = file_contents[seq(line_number+6,line_number+6+number_of_blocks)],strip.white=TRUE)
  colnames(raw_data)<-column_names
  
  isotopes <- column_names[3:length(column_names)]
  
  # Parse Primart Current START
  line_number <- grep(file_contents,pattern="Primary Current START ")
  primary_current_start<-str_match(file_contents[line_number],"[:digit:].*")
  primary_current_start <- as.numeric(primary_current_start)
  
  
  # Parse isotopic ratios
  first_line_number <- grep(file_contents,pattern="ISOTOPICS") + 1
  last_line_number <- grep(file_contents,pattern="STATISTICS RESULTS") - 1
  isotopic_ratios<-file_contents[first_line_number:last_line_number];
  pos<-str_locate(isotopic_ratios,"\\[.*\\[")[,2] # Find lines where we have ratios of isotopes
  isotopic_ratios<-isotopic_ratios[which(pos>0)]
  temp<-as.data.frame(str_split(isotopic_ratios,"\\/",simplify=TRUE)); # split either side of the division sign
  left<-str_match(temp$V1,"\\*.*\\[");
  left<-str_sub(left,2,str_length(left)-2);
  right<-str_match(temp$V2,"\\*.*\\[");
  right<-str_sub(right,2,str_length(right)-2);
  ratio_expression<-paste0("raw_data$`",left,"`/raw_data$`",right,"`");
  
  for (i in 1:length(ratio_expression)) {
    evaluated_ratio_expression <- eval(parse(text=ratio_expression[i]));
    mean_evaluated_ratio <- mean(evaluated_ratio_expression);
    se_evaluated_ratio <- sd(evaluated_ratio_expression)/sqrt(length(evaluated_ratio_expression));
  }
  
  
  
  # combine desired data in a row of dataframe
  # this can then be combined with each successive dataframe as files are read in 
  # this is not the most efficient way to do this but it works, which is good.
   read_row <- data.frame(spotname, session_datetime, primary_current_start, mean_evaluated_ratio, se_evaluated_ratio, x_position, y_position)
  

  combined_df <- rbind(combined_df, read_row)
}
Sys.setenv(TZ='GMT')
V_SMOW <- 0.0020052
# properly name columns of dataframe and create sample column
colnames(combined_df) <- c("Spot", "DateTime", "Primary_current", "O18_O16", "O18_O16_1se", "X_pos", "Y_pos")
combined_df$sample <- as.factor(gsub('\\@.*','', sub('.*\\-', '', combined_df$Spot)))

combined_df <- combined_df %>% mutate(delta18 = ((O18_O16/V_SMOW)-1)*1000)
combined_df <- combined_df %>% mutate(delta18_1se = (O18_O16_1se/O18_O16)*1000)


## UI function definition

drift_corr_func <- function(combined_df){
  
  # User interface ----
  ui <- miniPage(
    gadgetTitleBar("Drift correction"),
    
    miniTabstripPanel(
      miniTabPanel("Data filtering", icon = icon("map-o"),
                   miniContentPanel(
                     fillCol(flex = c(2,1,2),
                          plotOutput("O_ratio_plot", height = "100%", click = "point_click", brush = "points_brushed"),
                          fillRow(flex = c(1,2),
                            fillCol(
                            checkboxInput("drift_corr_check", label = "Drift correction", value = TRUE),
                            textOutput(outputId = "drift_coefficients")),
                            sliderInput("Date_range_selector", "Select Date Range",
                                        min = min(as.POSIXct(combined_df$DateTime)),
                                        max = max(as.POSIXct(combined_df$DateTime)),
                                        value = c(min(as.POSIXct(combined_df$DateTime)), max(as.POSIXct(combined_df$DateTime))))
                          ),
                          plotOutput("Ip_plot", height = "100%")
                     #tableOutput("values"),

                   )),        
                   miniButtonBlock(
                     actionButton("add", "", icon = icon("thumbs-up")),
                     actionButton("sub", "", icon = icon("thumbs-down")),
                     #actionButton("none", "", icon = icon("ban")),
                     actionButton("all", "", icon = icon("refresh"))
                   )
      )
    )
  )
  
  
  # Server logic ----
  server <- function(input, output) {
    
    vals <- reactiveValues(keep = rep(TRUE, nrow(combined_df))) 
    
    
    
output$O_ratio_plot <- renderPlot({
      
        p <- ggplot(combined_df %>% filter("91500" == sample), aes(x = DateTime, y = delta18)) +
        geom_point() +
        geom_point(data = combined_df[!vals$keep, , drop = FALSE], colour = "grey80") +
        geom_errorbar(aes(ymin = delta18-(2*delta18_1se), ymax = delta18+(2*delta18_1se)))
      
      if(input$drift_corr_check == TRUE){
        p +
        geom_smooth(data = combined_df[vals$keep, , drop = FALSE] %>% filter("91500" == sample) %>% filter_by_time(.date_var = DateTime, .start_date = as.POSIXct(input$Date_range_selector[1]), .end_date = as.POSIXct(input$Date_range_selector[2])), formula = y ~ x, method = lm)
      }
      
      else if(input$drift_corr_check == FALSE){
      p
      }
        
      else(NULL)
    })
    

    
    
drift_coefficientslm <- reactive({
    lm(delta18~as.POSIXct(DateTime), data = combined_df[vals$keep, , drop = FALSE] %>% filter("91500" == sample) %>% filter_by_time(.date_var = DateTime, .start_date = as.POSIXct(input$Date_range_selector[1]), .end_date = as.POSIXct(input$Date_range_selector[2]))) 
    })

output$drift_coefficients <- renderText({
  if(input$drift_corr_check == TRUE){
  paste0("R^2 = ", round(summary(drift_coefficientslm())$adj.r.squared, digits = 3))
  }
  
  else if(input$drift_corr_check == FALSE){
  "No drift correction"
  }
  
  else(NULL)
  })


sliderValues <- reactive({
  
  data.frame(
    Name = c("Start","End"),
    Value = as.character(c(as.POSIXct(input$Date_range_selector[1]),as.POSIXct(input$Date_range_selector[2]))),
    stringsAsFactors = FALSE)
  
})

# Show the values in an HTML table ----
output$values <- renderTable({
  sliderValues()
})
    output$Ip_plot <- renderPlot({
      
        ggplot(combined_df, aes(x = DateTime, y = Primary_current)) +
        geom_point() +
        geom_point(data = combined_df %>% filter_by_time(.date_var = DateTime, .start_date = as.POSIXct(input$Date_range_selector[1]), .end_date = as.POSIXct(input$Date_range_selector[2])), colour = "red") +
        geom_point(data = combined_df[!vals$keep, , drop = FALSE], colour = "grey80")
      

    })
    

    
    selected <- reactive({
      #nearPoints(combined_df, input$point_click, maxpoints = 1, allRows = TRUE)$selected_
      brushedPoints(combined_df, input$points_brushed, allRows = TRUE)$selected_
    })
    
    
    observeEvent(input$add, vals$keep <- vals$keep | selected())
    observeEvent(input$sub, vals$keep <- vals$keep & !selected())
    observeEvent(input$all, vals$keep <- rep(TRUE, nrow(combined_df)))
    observeEvent(input$none, vals$keep <- rep(FALSE, nrow(combined_df)))
    
    observeEvent(input$done, {
      stopApp(vals$keep)
    })
  }
  
  
  # Run app ----
  runGadget(ui, server)
  
}


# run UI
drift_corr_func(combined_df)

drift_coefficientslm <- lm(delta18~as.POSIXct(DateTime), data = combined_df %>% filter("91500" == sample) %>% filter_by_time(.date_var = DateTime, .start_date = as.POSIXct("2021-02-12 01:53:00"), .end_date = as.POSIXct("2021-02-12 7:18:00"))) 
summary(drift_coefficientslm)


signif(summary(drift_coefficientslm)$adj.r.squared, 3)
summary(drift_coefficientslm)$coefficients[,1]


