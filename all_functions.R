#set working directory
setwd('C:/Users/Avery/OneDrive/Documents/2022 Internship')
#load libraries
pacman::p_load(GEOquery, limma, umap, ggplot2, dplyr, tibble, readr, ggrepel)

## isolate columns and change column names for GEO full table ##
get.results <- function(gse.data){
  log2FoldChange <- gse.data$logFC
  ID <- gse.data$ID
  pvalue <- gse.data$P.Value
  padj <- gse.data$adj.P.Val
  rownum <- row.names(gse.data)
  results = cbind(ID, log2FoldChange, pvalue, padj, rownum)
  results = as.data.frame(results)
}

## change all numbers from character to numeric ##
numeric <- function(results){
  columns <- c('log2FoldChange', 'pvalue', 'padj', 'rownum')
  results[, columns] <- lapply(columns, function(x) as.numeric(results[[x]]))
  return(results)
}

## identify results as diffexpressed or not using decideTests ## 
how.regulated <- function(results){
  #decideTests needs matrix of pvalues
  pvalues <- as.matrix(results$pvalue)
  #find which genes are significant w/p-value cut-off of .05
  sig <- decideTests(pvalues, adjust.method = "fdr", p.value=0.05) 
  #change to df for merging later
  sig <- as.data.frame(sig)
  colnames(sig)[1] <- "type"
  #need a column to merge by
  sig$rownum <- rownames(sig)
  #add IDs to data by merging w/reformatted GEO table
  ultimate <- merge(results, sig, by = "rownum")
  return(ultimate)
}

## separate up from down-regulated and color-code ##
color.code <- function(ultimate){
  ultimate$diffexpressed[ultimate$log2FoldChange>0 & ultimate$type == 1] <- "Up" #up-regulated
  ultimate$diffexpressed[ultimate$log2FoldChange<0 & ultimate$type == 1] <- "Down" #down-regulated
  #remove probe.names, type
  ultimate <- subset(ultimate, select = -c(type, rownum))
  #labels needed for ggplot
  ultimate$labels <- NA
  return(ultimate)
}

## return the list of the IDs for DEGs ##
return.DEGs <- function(ultimate) {
  list.DEGs <- ultimate[which(ultimate$diffexpressed == 'Up'|
                                ultimate$diffexpressed == 'Down'),]
  list.DEGs <- subset(list.DEGs, select = c('ID', 'diffexpressed', 'padj'))
  return(list.DEGs)
}

## make volcano plot ##
volcano.plot <- function(results, gse){
  volcano = ggplot(data = results, 
                   aes(x= log2FoldChange, 
                       y = -1*log10(pvalue), 
                       col = diffexpressed, 
                       label = labels))
  plot <- volcano+
    geom_point()+
    theme_minimal()+
    geom_text_repel()+
    #color-code up vs. down genes based on for.colors function
    scale_color_manual(values = c('blue','red'))+ 
    theme(text=element_text(size=14))+ #change text size
    ggtitle(gse)+ 
    xlab("log2(fold change)")+ 
    ylab("-log10(pvalue)")
  return(plot)
}

## convert IDs to Ensembl ##
require('biomaRt')
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'affy_hg_u133_plus_2',
    'ensembl_gene_id'),
  uniqueRows = TRUE)
names(annotLookup)[names(annotLookup) == 'affy_hg_u133_plus_2'] <- "ID"

#change names of diffexpressed and padj, then get rid of og label
names(T3.conv)[names(T3.conv) == 'diffexpressed'] <- "T3.diff"
names(T3.conv)[names(T3.conv) == 'padj'] <- "T3.padj"
T3.conv <- subset(T3.conv, select = c('ensembl_gene_id', 'T3.diff', 'T3.padj'))
#merge unique list of CAD and one of T2D
all.c <- merge(C1.conv, C2.conv, by = 'ensembl_gene_id', all = TRUE)


## resolve conflicting regulation ##
#find min
library(matrixStats)
all.c$min = rowMins(as.matrix(all.c[, c("C1.padj", "C2.padj")]))
all.t$min <- apply(all.t[, c("T1.padj", "T2.padj", "T3.padj")], 1, min, na.rm = TRUE)
#if-else for CAD column
all.c$final <- with(all.c, ifelse(is.na(C1.diff), C2.diff,
                                  ifelse(is.na(C2.diff), C1.diff,
                                         ifelse(min == C1.padj, C1.diff,
                                                ifelse(min == C2.padj, C2.diff, NA)))))
#incomplete for T2D
all.t$final <- with(all.t, ifelse(min == T1.padj, T1.diff,
                                  ifelse(min == T3.padj, T3.diff, T2.diff)))

#subset only final regulation and ID
c.final <- subset(all.c, select = c('ensembl_gene_id', 'final'))