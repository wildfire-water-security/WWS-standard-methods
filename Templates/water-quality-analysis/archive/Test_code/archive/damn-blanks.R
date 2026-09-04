#optional: add grab samples. Not sure how to code this yet.How to include in random numbers???? Make own dataframe? Just have blank IDs in the lab?



GRAB <- FALSE #fill TRUE for having grab samples, FALSE for not. 
grabsamp <- c(2,0) #grab samples per site 

pbn <- 4 #number of procedural blanks 
fb <- 
  
  
  blanks <- (procedural_blank = "procedural_blank")
blanks <- rep(blanks, times = procedural_blank)
make.unique(blanks, sep = "_")

extra <- data.frame(
  study =study,
  site= rep(sampleinfo$site, times = grabsamp),
  datetime_PST= blanktime,
  samp_type = "GB"
)
problank

"PB"

#will need to only filter out IS for field labels??????????












blanks
combined_df <- bind_rows(data, blanks)
randocom <- combined_df %>%
  mutate (randomn= sample(rn:ts,ts, replace=F))
randocom

bs <- sum(grabsamp)
ts <- sum(sn,bs
          
)
ts
#extras
extras <- TRUE

field_blanks <- TRUE
grab_samples <- TRUE

grabsites <- c(STA,COA)
grabnumber <- c()



pb <- 4 #suggested to do one for each campaign and extra for each 50-100?
pro_BLK <- "pro_BLK" # naming for procedural blanks

num = (1:pb)
blk <- data.frame(
  samp_type = "PB",
  study=study,
  name = paste0(pro_BLK,"_",num),
  blanktime =blanktime
)
blk2 <-blk %>%
  mutate(
    field_sample_ID=paste0(study, "_",
                           name, "_",
                           samp_type,"_",
                           blanktime))
blk2

#SSEMG_pro_BLK_1

ts <- sum(sn,pb)
ts
combined_df <- bind_rows(data, blk)

combined_df



randocom <- combined_df %>%
  mutate (randomn= sample(rn:ts,ts, replace=F))
randocom

if(combined_df$samp_type == "PB" ){
  datalabels3<-combined_df %>%
    mutate(field_sample_ID = paste0(study, "_",
                                    site, "_",
                                    samp_type,"_",
                                    datetime_PST))
}



