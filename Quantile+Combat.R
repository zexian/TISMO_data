library(optparse)
library(tidyr)
library(reshape2)
library(dplyr)
library(preprocessCore)
library(limma)


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


data <- readRDS(opts$input)
meta <- read.csv(opts$meta, stringsAsFactors = FALSE)

cline <- unique(as.character(meta$Cell_Line))
exp_quantile <- data[1]
colnames(exp_quantile) <- "tmp"

for (c in cline) {
  stu_num <- length(unique(meta[meta$Cell_Line == c,]$GSE_ID))
  sample_id <- as.character(meta[meta$Cell_Line == c,]$SampleName)
  tmp_exp <- data[colnames(data) %in% sample_id]
  print(paste0("processing the ", sample_id, sep = ""))
  if (stu_num > 1) {
    tmp_exp <- as.matrix(tmp_exp)
    tmp_quantiles <- normalize.quantiles(tmp_exp)
    rownames(tmp_quantiles) <- rownames(tmp_exp)
    colnames(tmp_quantiles) <- colnames(tmp_exp)
    tmp_quantiles <- as.data.frame(tmp_quantiles)
  } else {
    tmp_quantiles <- tmp_exp
  }
  exp_quantile <- cbind(exp_quantile, tmp_quantiles)
}

exp_quantile <- exp_quantile[-1]
#exp_loess_logform <- log2(exp_loess + 1)

#saveRDS(exp_quantile,opts$output)
print("Running quantile normalization successfully!")
print("Start applying combat")

#combat
expr.dat <- log2(exp_quantile)
###filtering out genes with low variance among samples
CVFILTER <- 0
mean_nolym <- apply(expr.dat,1,mean)
var_nolym <- apply(expr.dat,1,var)
cv_nolym <- abs(var_nolym/mean_nolym)
filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)

## Select those genes that pass variance filtering
exprZero <- expr.dat
expr.dat <- expr.dat[rownames(expr.dat) %in% names(filt_genes),]
exprZero <- subset(exprZero, !(rownames(exprZero) %in% names(filt_genes)))

#batch removal using GSE ID within same cell line 
covars <- "GSE_ID"
covars <- append(covars, "tmp")

#create a tmp final data frame
expr.combat.all <- expr.dat[1]
colnames(expr.combat.all) <- "tmp"

for (cl in unique(meta$Cell_Line)) {
  
  #subset the samples based on cell line
  meta_sub <- meta[meta$Cell_Line == cl,]
  meta_sub$GSE_ID <- factor(meta_sub$GSE_ID, levels = unique(meta_sub$GSE_ID))
  
  sampleid <- as.character(meta_sub$SampleName)

  if(length(unique(meta_sub$GSE_ID)) == 1) {
    expr.combat <- expr.dat[,sampleid]
  } else {
    #add temporary column to fit for model.matrix function
    rownames(meta_sub) <- meta_sub$SampleName
    meta_sub$tmp <- 1
    print(paste0("cell line to be adjusted:", cl))
    
    meta_sub <- meta_sub[sampleid,covars]
    meta_sub <- rbind(data.frame(GSE_ID = "tmp", tmp = 1), meta_sub)
    #add temporary column to fit for model.matrix function
    fomul <- as.formula(paste("~ ", paste0(covars, collapse = " + "), sep = " "))
    design <- model.matrix(fomul, data=meta_sub)
    
    #remove the tmp column
    design <- design[-1,]
    design <- design[,c(-1, -ncol(design))]
    
    #running combat
    expr.combat = tryCatch(
      removeBatchEffect(as.matrix(expr.dat[,sampleid]),
                        covariates = design),
      error = function(e){
        print(e)
      })
    expr.combat <- as.data.frame(expr.combat)
  }
 
  expr.combat.all <- cbind(expr.combat.all, expr.combat)
}

expr.combat.all <- expr.combat.all[-1]

expr.combat.all <- rbind(expr.combat.all, exprZero)

saveRDS(expr.combat.all,opts$output)
print("Finished!")

