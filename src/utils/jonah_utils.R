

# Helps remove outliers
clamp <- function(lst, max = Inf, min = -max, rm = FALSE) {
  if (rm == TRUE) {
    lst = lst[lst < max]
    lst = lst[lst > min]
    return(lst)
  }
  else {
    return(pmax(pmin(lst, max), min))
  }
}

# Faster import/export
mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}
mcreadRDS <- function(file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -d -c -p",mc.cores," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}