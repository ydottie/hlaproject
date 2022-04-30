

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

cal_acc <- function(dat, gsfile, type=4){
# cal_acc <- function(dat, type=2){
    
  gs <- fread(gsfile) %>% as_tibble
  gs <- gs %>% arrange(dataset)
  
  stopifnot(all.equal(pull(gs[,1]), pull(dat[,1])))
  
  # dat <- apply(dat, 2, function(x){
  #   sub("([A-C]*[0-9]{2}:[0-9]{2}(|[A-z])).*", "\\1", x)
  # }) %>% as_tibble
  
  if(type==2){
    for(i in 2:7){
      gs[,i]<- sub("([A-z]*[0-9]{2}).*","\\1",pull(gs[,i]))
    }
  }else{
    idx_xx <- data.frame(which(apply(gs, c(1,2), function(x) grepl(":xx", x)), arr.ind = TRUE))
    for(i in 1:length(idx_xx)){
      gs[idx_xx$row[i], idx_xx$col[i]] <- sub("([A-z]*[0-9]{2}).*","\\1",gs[idx_xx$row[i], idx_xx$col[i]])
      # dat[idx_xx$row[i], idx_xx$col[i]] <- sub("([A-z]*[0-9]{2}).*","\\1",dat[idx_xx$row[i], idx_xx$col[i]])
    }
  }
  
  idx_type <- list(2:3, 4:5, 6:7)
  names(idx_type) <- c("A","B","C")
  
  for(type in names(idx_type)){
    ABC_gs <- gs[,idx_type[[type]]]
    ABC_dat <- dat[,idx_type[[type]]]
    
    tmp_res <- sapply(1:nrow(ABC_gs), function(i){
      
      sel_ABC_gs <- unique(unlist(ABC_gs[i,]))
      sel_ABC_dat <- unique(unlist(ABC_dat[i,]))
      
      if(length(sel_ABC_gs) == 2 && length(sel_ABC_dat) == 2){
        sum(grepl(sel_ABC_gs[1], sel_ABC_dat, fixed = T), 
            grepl(sel_ABC_gs[2], sel_ABC_dat, fixed = T))
      }else if(length(sel_ABC_gs) == 2 && length(sel_ABC_dat) == 1){
        sum(grepl(sel_ABC_gs[1], sel_ABC_dat, fixed = T), 
            grepl(sel_ABC_gs[2], sel_ABC_dat, fixed = T))
      }else if(length(sel_ABC_gs) == 1 && length(sel_ABC_dat) == 2){
        sum(grepl(sel_ABC_gs[1], sel_ABC_dat, fixed = T))
      }else{
        sum(grepl(sel_ABC_gs, sel_ABC_dat, fixed = T)) * 2
        # 2
      }
    })
    
    if(type == "A") {
      res <- data.frame(tmp_res) 
    }else{
      res <- cbind(res, data.frame(tmp_res))  
    }
  }
  
  colnames(res) <- names(idx_type)
  
  for(x in 1:3){
    res[,x][res[,x] > 2] <- 2
  }
  
  final_res <- res %>% as_tibble %>% 
    mutate(fourdigit=rowSums(res)/6) %>%
    mutate(gs[,1]) %>% 
    select(dataset,everything())
  
  cat("Acuracy:", mean(final_res$fourdigit),"\n")
  return(final_res)
}