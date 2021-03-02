## Version 0.1 Stable isotope data reduction ##



file_pathname <- file.choose()
file_contents = readLines(file_pathname)


# Parse acquisition filename
line_number <- grep(file_contents,pattern="ACQUISITION FILE NAME")
position <- str_locate(file_contents[line_number],"\\t")
acquisition_filename <-  str_sub(file_contents[line_number],position[1]+1)
spotname <- sub("^.+\\\\", "", acquisition_filename)

# Parse session time
line_number <- grep(file_contents,pattern="CMCA")
session_time <- str_match(file_contents[line_number],"[:digit:].*")

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



# combine desired data in a vector
# this can then be combined with other vectors to form a table 
# each vector as a row
testrow <- c(spotname, session_time, primary_current_start, mean_evaluated_ratio, se_evaluated_ratio, x_position, y_position)


