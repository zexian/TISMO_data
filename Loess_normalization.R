library(optparse)
library(tidyr)
library(reshape2)

option_list <- list(
  make_option(c("-i", "--input"), type="character",
              help="path of exp RDS file. [Required]"),
  make_option(c("-m", "--meta"), type="character",
              help="path of meta. [Required]"),
  make_option(c("-o", "--output"), type="character",
              help="Output file names. [Required]")
)

opt_parser <- OptionParser(option_list=option_list,add_help_option = FALSE);
opts <- parse_args(opt_parser);

#loess normalization
normalize.loess <- function(mat, subset=sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
                            epsilon=10^-2, maxit=1, log.it=FALSE, verbose=TRUE, span=2/3,
                            family.loess="symmetric", ...){
  
  J <- dim(mat)[2]
  II <- dim(mat)[1]
  if(log.it){
    mat <- log2(mat + 1)
  }
  
  change <- epsilon +1
  iter <- 0
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
  ##extremes added so that we can interpolate
  
  while(iter < maxit){
    iter <- iter + 1
    means <- matrix(0,II,J) ##contains temp of what we substract
    
    for (j in 1:(J-1)){
      for (k in (j+1):J){
        y <- mat[,j] - mat[,k]
        x <- (mat[,j] + mat[,k]) / 2
        index <- c(order(x)[1], subset, order(-x)[1])
        ##put endpoints in so we can interpolate
        xx <- x[index]
        yy <- y[index]
        aux <-loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict(aux, data.frame(xx=x)) / J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
        if (verbose)
          cat("Done with",j,"vs",k,"in iteration",iter,"\n")
      }
    }
    mat <- mat - means
    change <- max(colMeans((means[subset,])^2))
    
    if(verbose)
      cat(iter, change,"\n")
    
  }
  
  if ((change > epsilon) & (maxit > 1))
    warning(paste("No convergence after", maxit, "iterations.\n"))
  
  if(log.it) {
    return(2^mat)
  } else
    return(mat)
}


data <- readRDS(opts$input)
meta <- read.csv(opts$meta, stringsAsFactors = FALSE)

cline <- unique(as.character(meta$Cell_Line))
exp_loess <- data[1]
colnames(exp_loess) <- "tmp"

for (c in cline) {
  stu_num <- length(unique(meta[meta$Cell_Line == c,]$GSE_ID))
  sample_id <- as.character(meta[meta$Cell_Line == c,]$SampleName)
  tmp_exp <- data[colnames(data) %in% sample_id]
  
  if (stu_num > 1) {
    tmp_loess <- normalize.loess(tmp_exp)
  } else {
    tmp_loess <- tmp_exp
  }
  exp_loess <- cbind(exp_loess, tmp_loess)
}

exp_loess <- exp_loess[-1]
#exp_loess_logform <- log2(exp_loess + 1)

saveRDS(exp_loess,opts$output)
print("Running successfully!")

