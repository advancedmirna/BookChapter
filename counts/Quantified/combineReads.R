###############################
# Author: Hrishikesh Lokhande
# Date: 11/30/2021
# Description: This program is provided with the book chapter
# entitied 'Adavanced Bioinformatics Analysis of MiRNA seq data'
# This program combines provides two different functions
# 1) Add same mirna precursors counts (Focus is only precursor DE)
# 2) Process individual read count file and produce a single
# expression matrix with all samples reads as columns.
# This program is run before running the DE analysis with DESeq2
#############################

library(dplyr)
unique_precursor <-
  read.table('total_precursor.txt', header = T) ## A unique file with all miRNA precursor names

files <-
  list.files(
    path = ".",
    pattern = "*.tsv",
    full.names = TRUE,
    recursive = FALSE
  )   ## Reads all files with paths with *.tsv extension into a list


expressionFrame <- lapply(files, function(x) {
  sample_name <- sub("*./", "", x)
  sample_name <-
    sub('.tsv', '', sample_name)  ## This ensures sample names are preserved and added later
  
  mir_frame <-
    read.csv(x, sep = "\t", header = TRUE) ## Reading individual files
  mir_frame <- mir_frame %>% select(precursor, read_count) %>%
    group_by(precursor) %>%
    summarize(x = sum(read_count)) %>%
    select(x)
  colnames(mir_frame) <- c(sample_name)
  return(mir_frame)
}) %>% bind_cols()  ## Binding of columns into a single dataframe

expressionFrame <- data.frame(expressionFrame)
rownames(expressionFrame) <- unique_precursor$Precursor


expressionFrame <-
  expressionFrame[rowSums(expressionFrame[]) > 0,] ## removing unexpressed miRNA precursors.

write.csv(expressionFrame, file = 'ExpressionMatrix.csv')
