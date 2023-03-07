aclnrintp <- function(x, threshold = 0.8){
    autocorr <- (acf(x[!is.na(x)]))$acf
    cut <- rle(as.vector(autocorr > threshold))
    cutoff <- cut$lengths[cut$values][1]
    nalen <- rle(as.vector(is.na(x)))
    crx <- x
    id <- 0
    for(i in 1:length(nalen$values)){
        if(nalen$values[i] == FALSE){
            id <- id + nalen$lengths[i]
        }else{
            if((nalen$lengths[i] < cutoff) & (id > 0)){
                lv <- x[id]
                hv <- x[id+nalen$lengths[i]+1]
                crx[(id+1):(id+nalen$lengths[i])] <- lv+((hv-lv)*(1:nalen$lengths[i])/(nalen$lengths[i]+1))
            }
            id <- id + nalen$lengths[i]
        }
    }
    return(crx)
}
