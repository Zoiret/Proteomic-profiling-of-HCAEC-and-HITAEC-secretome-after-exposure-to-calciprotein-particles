rm(list = ls())

{
library(BiocManager)
library(readxl)
library(devtools)
library(ggvenn)
library(impute)
library(RColorBrewer)
library(limma)
library(vsn)
library(NMF)
library(ggplot2)
library(mixOmics)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(openxlsx)
}


#Uploading and preparing data

fact <- data.frame(read_excel("sample_info_131-135_vs_142-147.xlsx"))
  rownames(fact) <- fact[,1]
  fact <- fact[,-1]
  fact$Differentiation <- factor(fact$Differentiation, levels = c("HIT", "HCA"))
  fact$Differentiation
  table(fact$Differentiation)
  fact$Donor <- as.factor(fact$Donor)
  fact$Donor
  
  fact$group <- substr(fact$Original.index, 1, 1) #workmane of group incubation conditions
  fact$group <- as.factor(fact$group)
  
  fact$group_2 <- paste(fact$Differentiation, fact$group, sep = "_") #workname in registry without number
  fact$group_2 <- as.factor(fact$group_2)
  
  fact$Original.index <- paste(fact$Differentiation, fact$Original.index, sep = "_") #workname in registry
  fact$Original.index <- as.factor(fact$Original.index)
  
  rownames(fact) <- fact$Original.index
  str(fact)

dat <- data.frame(read.delim("proteins_131-153.tsv"))   
  dat[dat==0] <- NA 

  dat1 <- dat[,c(1,109:112,119:124)]
  head(dat1)
  str(dat1)
  
  rownames(dat1) <- dat1[,1]
  dat1 <- dat1[,-1]
  head(dat1)
  summary(dat1)
  
  colnames(dat1) <- c("HCA_K2", "HCA_K4", "HCA_M2", "HCA_M6", "HIT_K2", "HIT_K4", "HIT_K6", "HIT_M2", "HIT_M4", "HIT_M6")
  dat1 <- dat1[,c("HIT_K2", "HIT_K4", "HIT_K6", "HIT_M2", "HIT_M4", "HIT_M6", "HCA_K2", "HCA_K4", "HCA_M2", "HCA_M6")]
  summary(dat1)

#Qualititative analysis

    #detection of optimal way is to filter out proteins with too much NA
    {
    df <- data.frame(part=0:1000, HCA_K=0, HCA_M=0, HIT_K=0, HIT_M=0)
    c <- matrix(,ncol = 1001, nrow=4)
    for (i in 0:1000) 
      {
      K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_K"))])) >= i/1000), ]
      M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_M"))])) >= i/1000), ]
      KK <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_K"))])) >= i/1000), ]
      MM <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_M"))])) >= i/1000), ]
      c[1,i+1] <- nrow(K)
      c[2,i+1] <- nrow(M)
      c[3,i+1] <- nrow(KK)
      c[4,i+1] <- nrow(MM)
      }
    df['HCA_K'] <- c[1,]
    df['HCA_M'] <- c[2,]
    df['HIT_K'] <- c[3,]
    df['HIT_M'] <- c[4,]
    ggplot(data = df) +
      #geom_smooth(mapping = aes(x = part, y = HCA_K, colour = "HCA_K")) + 
      geom_line(mapping = aes(x = part, y = HCA_K, colour = "HCA_K")) +
      #geom_smooth(mapping = aes(x = part, y = HCA_M, colour = "HCA_M")) + 
      geom_line(mapping = aes(x = part, y = HCA_M, colour = "HCA_M")) +
      #geom_smooth(mapping = aes(x = part, y = HIT_K, colour = "HIT_K")) + 
      geom_line(mapping = aes(x = part, y = HIT_K, colour = "HIT_K")) +
      #geom_smooth(mapping = aes(x = part, y = HIT_M, colour = "HIT_M")) + 
      geom_line(mapping = aes(x = part, y = HIT_M, colour = "HIT_M")) +
      geom_vline(xintercept = 600) +
      xlab("Decline NA (‰)") + ylab("Detected proteins") +
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      scale_y_continuous(breaks = seq(0, 3000, 500)) +
      labs(subtitle="Dependence of the number of detectable proteins on decline NA (HCA ctrl vs HIT ctrl)")
    }
 

    {
      df <- data.frame(part=0:1000, HCA=0, HIT=0)
      c <- matrix(,ncol = 1001, nrow=2)
      for (i in 0:1000) 
      {
        HCA <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HCA"))])) >= i/1000), ]
        HIT <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HIT"))])) >= i/1000), ]
        c[1,i+1] <- nrow(HCA)
        c[2,i+1] <- nrow(HIT)
      }
      df['HCA'] <- c[1,]
      df['HIT'] <- c[2,]
      ggplot(data = df) +
        geom_smooth(mapping = aes(x = part, y = HCA, colour = "HCA")) + 
        geom_line(mapping = aes(x = part, y = HCA, colour = "HCA")) +
        geom_smooth(mapping = aes(x = part, y = HIT, colour = "HIT")) + 
        geom_line(mapping = aes(x = part, y = HIT, colour = "HIT")) +
        geom_vline(xintercept = 600) +
        xlab("Decline NA (‰)") + ylab("Detected proteins") +
        scale_x_continuous(breaks = seq(0, 1000, 100)) +
        scale_y_continuous(breaks = seq(0, 3000, 500)) +
        labs(subtitle="Dependence of the number of detectable proteins on decline NA (HCA ctrl vs HIT ctrl)")
    }
 
  #filter out proteins with too much NA
  HCA_K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_K"))])) >= 0.5), ]
  HCA_M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_M"))])) >= 0.5), ]
  HIT_K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_K"))])) >= 0.5), ]
  HIT_M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_K"))])) >= 0.5), ]
  
  HCA <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HCA"))])) >= 0.5), ]
  HIT <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HIT"))])) >= 0.5), ]
  
  vennn <- list(HCA_K = rownames(HCA_K), HCA_M = rownames(HCA_M), HIT_K = rownames(HIT_K), HIT_M = rownames(HIT_M))
ggvenn(vennn, 
       fill_color = c("#0073C2FF", "#CD534CFF", "#009E73", "#F0E442"),
       stroke_size = 0.5, set_name_size = 8, text_size = 5,)
  
  vennn2 <- list(HCA = rownames(HCA), HIT = rownames(HIT))
ggvenn(vennn2, 
       fill_color = c("#0073C2FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 8, text_size = 5,)

  vennn3 <- list(HCA_K = rownames(HCA_K), HIT_K = rownames(HIT_K), HIT_M = rownames(HIT_M))
ggvenn(vennn3, 
       fill_color = c("#0073C2FF", "#CD534CFF", "#009E73", "#F0E442"),
       stroke_size = 0.5, set_name_size = 8, text_size = 5,)

  # Multiple set version of intersect
  Intersect <- function (x) {
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
      intersect(x[[1]], Intersect(x[-1]))
    }
  }
  
  # Multiple set version of union
  Union <- function (x) {  
    # x is a list
    if (length(x) == 1) {
      unlist(x)
    } else if (length(x) == 2) {
      union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
      union(x[[1]], Union(x[-1]))
    }
  }
  
  # Remove the union of the y's from the common x's
  Setdiff <- function (x, y) {
    # x and y are lists of characters
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
  }
  
  HCA_K_spec <- Setdiff(vennn[c("HCA_K")], vennn[c("HIT_M", 'HCA_M', 'HIT_K')])
  HCA_M_spec <- Setdiff(vennn[c("HCA_M")], vennn[c("HIT_M", 'HCA_K', 'HIT_K')])
  HIT_K_spec <- Setdiff(vennn[c("HIT_K")], vennn[c("HIT_M", 'HCA_M', 'HCA_K')])
  HIT_M_spec <- Setdiff(vennn[c("HIT_M")], vennn[c("HCA_K", 'HCA_M', 'HIT_K')])
  
  HCA_spec <- Setdiff(vennn2[c("HCA")], vennn2[c("HIT")])
  HIT_spec <- Setdiff(vennn2[c("HIT")], vennn2[c("HCA")])

#   write.table(HCA_K_spec, "HCA_K_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
#   write.table(HCA_M_spec, "HCA_M_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
#   write.table(HIT_K_spec, "HIT_K_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
#   write.table(HIT_M_spec, "HIT_M_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
#   write.table(HCA_spec, "HCA_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
#   write.table(HIT_spec, "HIT_spec_131-135_vs_142-147.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

#Quantitative analysis

  dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.71),] #filtration data according to the above algorithm
  mean(complete.cases(dat2))
  NAsums <- data.frame(colSums(is.na(dat2)))
  NAsums
  # write.table(NAsums, "NA_sums_142-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  str(dat2)

#Imputation
  tdat <- t(dat2)
  dat_knn1 <- impute.knn(tdat, k = 5)
  dat_knn <- t(dat_knn1$data)
  mean(complete.cases(dat_knn))

#Structure of row data
  pal <- brewer.pal(n = 9, name = "Set1")
  cols <- pal[fact$group_2]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data (HCA ctrl vs HIT ctrl)")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)
  colSums(dat_knn)
  
  length(dat_knn[dat_knn == 0]) 

#Logarithm of data
  dat_log <- log2(dat_knn+1) 
  head(dat_log)
  mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data (HCA ctrl vs HIT ctrl)")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)

#Normalization of ata
  dat_norm <- normalizeQuantiles(dat_log) 
  head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data (HCA ctrl vs HIT ctrl)")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)
  mean(complete.cases(dat_norm))
  colSums(is.na(dat_norm))

#MAplot (Log-expression)
  maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
    X <- (rowMeans(X2) + rowMeans(X1)) / 2
    Y <- rowMeans(X2) - rowMeans(X1)
    scatter.smooth(x = X, y = Y,
                   main = main, pch = pch,
                   xlab = xlab, ylab = ylab,
                   lpars = lpars, ...)
    abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
  }

maplot(dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], 
       dat_log[, rownames(fact)[fact$Differentiation == "HIT"]], 
       main = "Log-expression data HCA ctrl vs HIT ctrl")
maplot(dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], 
       dat_norm[, rownames(fact)[fact$Differentiation == "HIT"]], 
       main = "Normalized data HCA ctrl vs HIT ctrl")

#MeanSd
meanSdPlot(as.matrix(dat_log))
meanSdPlot(as.matrix(dat_norm))

#Heatmap
aheatmap(cor(dat_norm), color = "-RdBu:256", annCol = fact$group_2, fontsize = 10)

#Principle components analysis
  dat_pca <- pca(t(dat_norm), ncomp = 10, center = TRUE)
  dat_pca
  
plot(dat_pca)

plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$Differentiation, legend = TRUE, ellipse = T,
          title = 'PCA HCA ctrl vs HIT ctrl') 

plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$group_2, legend = TRUE, ellipse = T,
          title = 'PCA HCA ctrl vs HIT ctrl') 

#Limma - differentially expressed proteins

  X <- model.matrix(~ fact$Differentiation)
  X
  
  fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)
  
  efit <- eBayes(fit)
  
  topTable(efit, coef = 2)
  numGenes <- length(dat_norm)
  full_list_efit <- topTable(efit, number = length(dat_norm[,1]))
  # write.csv(full_list_efit,'Dif_expr_HCA_ctrl_vs_HIT_ctrl.csv')
  head(full_list_efit)
  
  full_list_efit$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit$Gene.name)) 
  full_list_efit$Gene.name[is.na(full_list_efit$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit)[is.na(full_list_efit$Gene.name)])

  full_list_efit['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit))
  full_list_efit['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit))
  full_list_efit <- full_list_efit[,c(7,8,9,1,2,4,5)]
  
    # write.xlsx(full_list_efit, file = "CTRLs_HCS_vs_HIT_diff_proteins.xlsx")
  
#Vulcano plot
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-3.5, 3.5), 
                ylim = c(0, 5.5),
                FCcutoff = 1,  
                title ="Vulcano HCA ctrl vs HIT ctrl",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)



#Heatmap

  #rename values for building modified heatmap legend
  clusters_col <- c('HIT_K2', 'HIT_K4', 'HIT_K6', 'HIT_M2', 'HIT_M4', 'HIT_M6', 'HCA_K2', 'HCA_K4', 'HCA_M2', 'HCA_M6')
  print(colnames(dat2))
  print(colnames(dat_norm))
  dat2_2 <- dat2[,clusters_col]
  dat_norm2 <- dat_norm[,clusters_col]
  print(colnames(dat2_2))
  print(colnames(dat_norm2))
  print(rownames(dat_norm2))
  
  full_list_efit2 <- full_list_efit
  full_list_efit3 <- full_list_efit[full_list_efit$'adj.P.Val'< 0.05,]
  full_list_efit3 <- full_list_efit3[order(full_list_efit3$logFC,decreasing = T),]
  full_list_efit3 <- full_list_efit3[- grep("ALB", full_list_efit3$Gene.name),]
  
  rownames(full_list_efit2) <- paste(full_list_efit$Name.of.protein, full_list_efit$Gene.name, sep = '--')
  rownames(full_list_efit3) <- paste(full_list_efit3$Name.of.protein, full_list_efit3$Gene.name, sep = '--')
  rownames(dat2_2) <- paste(full_list_efit2$Name.of.protein,full_list_efit2$Gene.name, sep = '--')
  rownames(dat_norm2) <- paste(full_list_efit2$Name.of.protein,full_list_efit2$Gene.name, sep = '--')
  
  #selected proteins for heatmap visualization. They be selected on the basis of average expression 
  p_above_FC_top_for_heatmap <- rownames(full_list_efit3[c('P62807','O00231','P15559','P02774','P10124','P10909','Q8WZ75','P02788','Q6YHK3','P07998', #top10 of overexpression
                                                          'Q8IV08','Q04837','O00469','P08603','P13686','P37840','P15090','P21810','Q13201','P69905'),]) #top10 of gipoexpression
  length(p_above_FC_top_for_heatmap)
  
pheatmap(dat_norm2[p_above_FC_top_for_heatmap,],
         main = "Heatmap HCA ctrl vs HIT ctrl",
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = T,
         cluster_rows = F,
         cellwidth = 10,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff",'#005aeb','#240935','#b300b3',"#e600e6"))(100))

  min(dat_norm2[p_above_FC_top_for_heatmap,])
  max(dat_norm2[p_above_FC_top_for_heatmap,])
  
  