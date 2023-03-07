rmNA <- function(val, dfoar){
  nc <- dim(dfoar)[2]
  if(is.null(dim(dtdzslpas)[1])){
    dfoar[dfoar == val] <- NA
  }else{
    for(i in 1:nc){
      dfoar[dfoar[,i] == val,i] <- NA
    }
  }
  return(dfoar)
}
