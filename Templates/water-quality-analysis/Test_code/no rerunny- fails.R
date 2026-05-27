


loaded <- FALSE  #if the environment is still loaded or not



savesetting <- data.frame(site=site,nsamp=nsamp,interval=interval,samp_type=samp_type,blanktime=blanktime,sc=sc,study=study,TZone=TZone,
                          DOM=DOM,DOC=DOC,CCAL=CCAL,Levoglucosan=Levoglucosan,Excess=Excess,
                          fileformat=fileformat,blank_fileformat=blank_fileformat) 

savesetting <- data.frame(site,nsamp,interval,samp_type,sc,study,TZone,sn,
                          DOM,DOC,CCAL,Levoglucosan,Excess,
                          fileformat,blank_fileformat) 
savesetting
write.csv(savesetting,"savesetting.csv")





#if we wanted to not have people rerun it but it too hard :c
if(loaded==FALSE){
  saved <- read.csv("savesetting.csv") #bring in the saved settings
  site <-c(saved$site)
  interval <- c(saved$interval)
  samp_type  <- saved$samp_type
  sc  <- saved$sc
  sn <- saved$sn
  study <- saved$study
  TZone <- saved$TZone
  DOM <- saved$DOM
  DOC <- saved$DOC
  CCAL <- saved$CCAL
  Levoglucosan  <- saved$Levoglucosan 
  Excess  <- saved$Excess 
  blank_fileformat <-saved$blank_fileformat
  fileformat <- saved$fileformat
}  
saved