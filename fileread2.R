##To read all the files of throughfall and GW in the folder
fileread2 <- function(path2dat){
  setwd(path2dat)
  file_list <- list.files(path = path2dat, pattern = "\\.CSV$")
  comlist <- list(NA)
  for(i in 1:length(file_list)){
      tempdataset <- read.csv(file_list[i])
          if(ncol(tempdataset) == 5){
              tempdataset$X <- NA
          }
      tempdataset$Name <- file_list[i]
      comlist[[i]] <- tempdataset
  }
  return(comlist)
}
