#Test for label maker code? 
#import datasheets 
labsheet <- read.csv("labsamples_SSEMG_IS_20260223.csv")
labsheet
samplesheet <- read.csv("fieldsamples_SSEMG_IS_20260223.csv")
samplesheet

library(dplyr)
library(stringr)
library(kableExtra)
#Filtering sheet: organize by random IDs, include filters and analyte bottles
#fieldsheet: organize by samples, include extra columns 



kable(samplesheet, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "hold_position"), 
                full_width = FALSE, 
                font_size = 10) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#D7261E")
