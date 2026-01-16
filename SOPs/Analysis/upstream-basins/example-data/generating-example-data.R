# script used to generate the example data used for running upsream-basins.rmd

#create sample .csv
  #used example EWEB sites
  df <- data.frame(siteID = c("E020", "E470", "E543"),
                   locality = c("McKenzie River", "Quartz Creek", "Lookout Creek"),
                   latitude = c(44.05611154,44.12330799, 44.21006528),
                   longitude = c(-122.8290005,-122.3786855,-122.2573696))
  write.csv(df, "SOPs/Analysis/upstream-basins/example-data/example-sites.csv", row.names = FALSE, quote = FALSE)
