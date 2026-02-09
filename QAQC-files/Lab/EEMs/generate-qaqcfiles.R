#used to create and save MDL and check standard files to a user's computer for the
 #qaqc checks for eemanalyzeR (https://katiewampler.github.io/eemanalyzeR/)

#DIRECTIONS:
  # 1. ensure you are in the WWS-standard-methods project
  # 2. run the following lines of code

library(eemanalyzeR)

#make mdl file for eemanalyzeR
create_mdl("QAQC-files/Lab/EEMs/long-term-standards/blanks",
           recursive = TRUE, type="eem", iblank="_blank")

create_mdl("QAQC-files/Lab/EEMs/long-term-standards/blanks",
           recursive = TRUE, type="abs", iblank="_blank")


#make std file for eemanalyzeR
create_std("QAQC-files/Lab/EEMs/long-term-standards/tea-standards", recursive = TRUE,
           type="eem", iblank="_blank", abs_pattern="Abs")

create_std("QAQC-files/Lab/EEMs/long-term-standards/tea-standards", recursive = TRUE,
           type="abs", abs_pattern="Abs")
