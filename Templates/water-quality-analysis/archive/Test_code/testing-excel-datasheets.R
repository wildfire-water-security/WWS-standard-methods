##testing how to do nice table formatting in excel
library(openxlsx)
library(openxlsx2)
library(dplyr)

#try to get styles from a template
ids <- wb_load("Templates/water-quality-analysis/Test_code/SSEMG_IS_20260223_label-table.xlsx") %>% wb_to_df(sheet = "lab_filter_sheets")
datasheets <- wb_load("Templates/water-quality-analysis/wqsample-label-datasheet-templates.xlsx")

datasheets <- datasheets %>% wb_add_data(sheet= "DS-Bottle Numbers", dims="A3", x=ids$sample_name)
datasheets <- datasheets %>% wb_add_data(sheet= "DS-ISCO Bottle Weight", dims="A3", x=ids$sample_name)
datasheets <- datasheets %>% wb_add_data(sheet= "DS-Filter Weight", dims="A3", x=ids$sample_name)
wb_save(datasheets, "Templates/water-quality-analysis/Test_code/datasheet-test.xlsx")
wb_open(datasheets)

sheet <- wb %>% wb_to_df(sheet="DS-Bottle Numbers")
sty <- wb_get_cell_style(wb, sheet="DS-Bottle Numbers", dims = "A3")

wb$add_worksheet("Sheet 1")

wb$set_cell_style(dims = "A1", sheet="Sheet 1", style=sty)

wb_open(wb)

