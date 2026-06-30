
###########################

library(readxl)
#for note sheets
labsheet <- read_excel("SSEMG_IS_20260224_label-table.xlsx", sheet ="filtering")
labsheet
labsheets <- data.frame(labsheet)
labsheets
#for weight
waterweight <- mutate(labsheets,pre_weight="", post_weight="",notes="")
waterweight
#for filtering
filtering <- mutate(labsheet, Aqualog = "",Shimadzu = "", CCAL ="",Levo="",Excess="",notes="")
filtering
filterweights <- mutate(labsheet,pre_weight="", post_weight="",notes="")
filterweights
#for labels in lab
lablabels <- read_excel("SSEMG_IS_20260224_label-table.xlsx",sheet = "lab_labels")
fieldlabels <-read_excel("SSEMG_IS_20260224_label-table.xlsx",sheet = "field_labels")

###########


#Test for label maker code? 

####bottles- fill out#####
DOM <- TRUE
DOC <- FALSE
CCAL <- TRUE
Levoglucosan <- TRUE
Excess <- TRUE
sn <- 5
#add analytes that are TRUE
analysis_opt <- c("Shimadzu", "Aqualog", "CCAL", "Excess", "Levoglucosan")
analysis_choice <- c(DOC, DOM, CCAL, Excess, Levoglucosan)
analysiscodes <- analysis_opt[analysis_choice]
analysiscodes 
#some kinda data frame for this 
analysisn <- length(analysiscodes) #how many analytes there are




data <- matrix(nrow=1,ncol=length(analysiscodes))
colnames(data) <- analysiscodes
data





#import datasheets 
samplesheet <- read.csv("example_label-table2.csv")
samplesheet

#how to add analytes automatically into a table?
analytes <- c(Aqualog)
analytes
labsheet <- mutate(samplesheet, Aqualog = "",Shimadzu = "", CCAL ="",Levo="",Excess="")
labsheet
library(dplyr)
library(stringr)
library(kableExtra)


#katie's example
kable(labsheet, booktabs = FALSE) %>%
  kable_styling(latex_options = c("striped", "hold_position","repeat_header"), 
                full_width = FALSE, 
                font_size = 10) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#D7261E")


samplesheet
 #GRAVEYARD
#add lines vertical
 # )%>%
  # column_spec(1:ncol(labsheet),border_left = T, border_right = T)
# left, center, center, right, right
# knitr::kable(iris2, align = "lccrr")


#https://github.com/EcologyR/labeleR/
#install.packages("labeleR")
#install.packages("tinytex")
library("labeleR")




