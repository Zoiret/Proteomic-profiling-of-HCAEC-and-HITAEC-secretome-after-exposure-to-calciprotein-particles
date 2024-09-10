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

fact <- data.frame(read_excel("sample_info_131-141.xlsx"))
  rownames(fact) <- fact[,1]
  fact <- fact[,-1]
  fact$Differentiation <- as.factor(fact$Differentiation)
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

  dat1 <- dat[,c(1,103:112)]
  head(dat1)
  str(dat1)

  rownames(dat1) <- dat1[,1]
  dat1 <- dat1[,-1]
  head(dat1)
  
  colnames(dat1) #after rename cols 'dat1' need verification of compliance names cols and rows
  rownames(fact)
  colnames(dat1) <- rownames(fact)


#Qualititative analysis

    #detection of optimal way is to filter out proteins with too much NA
    {
    df <- data.frame(part=0:1000, I=0, C=0, M=0, K=0)
    c <- matrix(,ncol = 1001, nrow=4)
    for (i in 0:1000) 
      {
      I <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="I"))])) >= i/1000), ]
      C <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="C"))])) >= i/1000), ]
      M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="M"))])) >= i/1000), ]
      K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="K"))])) >= i/1000), ]
      c[1,i+1] <- nrow(I)
      c[2,i+1] <- nrow(C)
      c[3,i+1] <- nrow(M)
      c[4,i+1] <- nrow(K)
      }
    df['I'] <- c[1,]
    df['C'] <- c[2,]
    df['M'] <- c[3,]
    df['K'] <- c[4,]
    ggplot(data = df) +
      # geom_smooth(mapping = aes(x = part, y = I, colour = "I")) + 
      geom_line(mapping = aes(x = part, y = I, colour = "I")) +
      # geom_smooth(mapping = aes(x = part, y = C, colour = "C")) + 
      geom_line(mapping = aes(x = part, y = C, colour = "C")) +
      # geom_smooth(mapping = aes(x = part, y = M, colour = "M")) + 
      geom_line(mapping = aes(x = part, y = M, colour = "M")) +
      # geom_smooth(mapping = aes(x = part, y = K, colour = "K")) + 
      geom_line(mapping = aes(x = part, y = K, colour = "K")) +
      geom_vline(xintercept = 600) +
      xlab("Decline NA (‰)") + ylab("Detected proteins") +
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      scale_y_continuous(breaks = seq(0, 3000, 500)) +
      labs(subtitle="Dependence of the number of detectable proteins on decline NA (HCA)")
    }

  #filter out proteins with too much NA
  C <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="C"))])) >= 0.5), ]
  K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="K"))])) >= 0.5), ]
  I <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="I"))])) >= 0.5), ]
  M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="M"))])) >= 0.5), ]
  
  vennn <- list(C = rownames(C), K = rownames(K), I = rownames(I))
ggvenn(vennn, 
       fill_color = c("#0073C2FF", "#CD534CFF", "#009E73"),
       stroke_size = 0.5, set_name_size = 8, text_size = 5,)

  vennn2 <- list(C = rownames(C), K = rownames(K), I = rownames(I), M = rownames(M))
ggvenn(vennn2, 
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
  
  
  C_spec <- Setdiff(vennn2[c("C")], vennn2[c("K", 'M', 'I')])
  M_spec <- Setdiff(vennn2[c("M")], vennn2[c("K", 'C', 'I')])
  I_spec <- Setdiff(vennn2[c("I")], vennn2[c("K", 'M', 'C')])
  K_spec <- Setdiff(vennn2[c("K")], vennn2[c("C", 'M', 'I')])
  # write.table(C_spec, "C_spec_131-141.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(I_spec, "I_spec_131-141.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(K_spec, "K_spec_131-141.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(M_spec, "M_spec_131-141.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

#Quantitative analysis

  dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.71),] #filtration data according to the above algorithm
  mean(complete.cases(dat2))
  NAsums <- data.frame(colSums(is.na(dat2)))
  NAsums
  #write.table(NAsums, "NA_sums_131-141.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  str(dat2)

#Imputation
  tdat <- t(dat2)
  dat_knn1 <- impute.knn(tdat, k = 5)
  dat_knn <- t(dat_knn1$data)
  mean(complete.cases(dat_knn))

#Structure of row data
  pal <- brewer.pal(n = 9, name = "Set1")
  cols <- pal[fact$group_2]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data (HCA)")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)
  colSums(dat_knn)
  
  length(dat_knn[dat_knn == 0]) 

#Logarithm of data
  dat_log <- log2(dat_knn+1) 
  head(dat_log)
  mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data (HCA)")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)

#Normalization of data
  dat_norm <- normalizeQuantiles(dat_log) 
  head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data (HCA)")
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

maplot(dat_log[, rownames(fact)[fact$group == "M"]], 
       dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Log-expression data HCA (HCA M)")
maplot(dat_norm[, rownames(fact)[fact$group == "M"]], 
       dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Normalized data HCA (HCA M)")

maplot(dat_log[, rownames(fact)[fact$group == "C"]], 
       dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Log-expression data HCA (HCA C)")
maplot(dat_norm[, rownames(fact)[fact$group == "C"]], 
       dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Normalized data HCA (HCA C)")

maplot(dat_log[, rownames(fact)[fact$group == "I"]], 
       dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Log-expression data HCA (HCA I)")
maplot(dat_norm[, rownames(fact)[fact$group == "I"]], 
       dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Normalized data HCA (HCA I)")

maplot(dat_log[, rownames(fact)[fact$group == "K"]], 
       dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Log-expression data HCA (HCA K)")
maplot(dat_norm[, rownames(fact)[fact$group == "K"]], 
       dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], 
       main = "Normalized data HCA (HCA K)")

#MeanSd
meanSdPlot(as.matrix(dat_log))
meanSdPlot(as.matrix(dat_norm))

#Heatmap
aheatmap(cor(dat_norm), color = "-RdBu:256", annCol = fact$group_2, fontsize = 10)

#Principle components analysis
  dat_pca <- pca(t(dat_norm), ncomp = 10, center = TRUE)
  dat_pca

plot(dat_pca)

plotIndiv(dat_pca, comp = c(1, 2), ind.names = FALSE, 
          group = fact$group, legend = TRUE, ellipse = TRUE,
          title = 'PCA HCA') 

#Limma - differentially expressed proteins
  design = model.matrix(~0+fact$group)
  colnames(design) <- c('C', 'I', 'K', 'M')
  colnames(design)
  
  fit <- lmFit(dat_norm, design=design, method = "robust", maxit = 10000)
  contrasts_groups = c('I-C','C-K','C-M','I-K','I-M','M-K') 
  contrast.matrix <- makeContrasts(contrasts = contrasts_groups,levels=design)
  
  efit <- contrasts.fit(fit, contrast.matrix)
  efit <- eBayes(efit)

  full_list_efit_IvsC <- topTable(efit, number = length(dat_norm[,1]), coef=1)
  full_list_efit_CvsK <- topTable(efit, number = length(dat_norm[,1]), coef=2)
  full_list_efit_CvsM <- topTable(efit, number = length(dat_norm[,1]), coef=3)
  full_list_efit_IvsK <- topTable(efit, number = length(dat_norm[,1]), coef=4)
  full_list_efit_IvsM <- topTable(efit, number = length(dat_norm[,1]), coef=5)
  full_list_efit_MvsK <- topTable(efit, number = length(dat_norm[,1]), coef=6)
  full_list_efit <- rbind(full_list_efit_IvsC, full_list_efit_CvsK, full_list_efit_CvsM, full_list_efit_IvsK, full_list_efit_IvsM, full_list_efit_MvsK)
  
  full_list_efit_IvsC$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_IvsC)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_IvsC$Gene.name)) 
  full_list_efit_IvsC$Gene.name[is.na(full_list_efit_IvsC$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_IvsC)[is.na(full_list_efit_IvsC$Gene.name)])

  full_list_efit_CvsK$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_CvsK)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_CvsK$Gene.name)) 
  full_list_efit_CvsK$Gene.name[is.na(full_list_efit_CvsK$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_CvsK)[is.na(full_list_efit_CvsK$Gene.name)])
  
  full_list_efit_CvsM$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_CvsM)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_CvsM$Gene.name)) 
  full_list_efit_CvsM$Gene.name[is.na(full_list_efit_CvsM$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_CvsM)[is.na(full_list_efit_CvsM$Gene.name)])
  
  full_list_efit_IvsK$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_IvsK)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_IvsK$Gene.name)) 
  full_list_efit_IvsK$Gene.name[is.na(full_list_efit_IvsK$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_IvsK)[is.na(full_list_efit_IvsK$Gene.name)])
  
  full_list_efit_IvsM$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_IvsM)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_IvsM$Gene.name)) 
  full_list_efit_IvsM$Gene.name[is.na(full_list_efit_IvsM$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_IvsM)[is.na(full_list_efit_IvsM$Gene.name)])
  
  full_list_efit_MvsK$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit_MvsK)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit_MvsK$Gene.name)) 
  full_list_efit_MvsK$Gene.name[is.na(full_list_efit_MvsK$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit_MvsK)[is.na(full_list_efit_MvsK$Gene.name)])
  
  full_list_efit$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit$Gene.name)) 
  full_list_efit$Gene.name[is.na(full_list_efit$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit)[is.na(full_list_efit$Gene.name)])
  
  full_list_efit_IvsC['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_IvsC))
  full_list_efit_CvsK['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_CvsK))
  full_list_efit_CvsM['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_CvsM))
  full_list_efit_IvsK['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_IvsK))
  full_list_efit_IvsM['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_IvsM))
  full_list_efit_MvsK['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit_MvsK))
  
  full_list_efit_IvsC['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_IvsC))
  full_list_efit_CvsK['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_CvsK))
  full_list_efit_CvsM['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_CvsM))
  full_list_efit_IvsK['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_IvsK))
  full_list_efit_IvsM['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_IvsM))
  full_list_efit_MvsK['Number.of.protein'] <- sub('.*\\|', '', rownames(full_list_efit_MvsK))
  
  full_list_efit_IvsC2 <- full_list_efit_IvsC[,c(7,8,9,1,2,4,5)]
  full_list_efit_CvsK2 <- full_list_efit_CvsK[,c(7,8,9,1,2,4,5)]
  full_list_efit_CvsM2 <- full_list_efit_CvsM[,c(7,8,9,1,2,4,5)]
  full_list_efit_IvsK2 <- full_list_efit_IvsK[,c(7,8,9,1,2,4,5)]
  full_list_efit_IvsM2 <- full_list_efit_IvsM[,c(7,8,9,1,2,4,5)]
  full_list_efit_MvsK2 <- full_list_efit_MvsK[,c(7,8,9,1,2,4,5)]
  
  
  # Differ_proteins <- createWorkbook()
  # 
  # addWorksheet(Differ_proteins, "HCAEC I-C")
  # writeData(Differ_proteins, "HCAEC I-C", full_list_efit_IvsC2)
  # 
  # addWorksheet(Differ_proteins, "HCAEC C-K")
  # writeData(Differ_proteins, "HCAEC C-K", full_list_efit_CvsK2)
  # 
  # addWorksheet(Differ_proteins, "HCAEC C-M")
  # writeData(Differ_proteins, "HCAEC C-M", full_list_efit_CvsM2)
  # 
  # addWorksheet(Differ_proteins, "HCAEC I-K")
  # writeData(Differ_proteins, "HCAEC I-K", full_list_efit_IvsK2)
  # 
  # addWorksheet(Differ_proteins, "HCAEC I-M")
  # writeData(Differ_proteins, "HCAEC I-M", full_list_efit_IvsM2)
  # 
  # addWorksheet(Differ_proteins, "HCAEC M-K")
  # writeData(Differ_proteins, "HCAEC M-K", full_list_efit_MvsK2)
  # 
  # saveWorkbook(Differ_proteins, file = "HCAEC_diff_proteins.xlsx", overwrite = TRUE)

#Vulcano plot
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano all vs all (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_IvsC,
                lab = rownames(full_list_efit_IvsC), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano I vs C (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_CvsK,
                lab = rownames(full_list_efit_CvsK), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano C vs K (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_CvsM,
                lab = rownames(full_list_efit_CvsM), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano C vs M (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_IvsK,
                lab = rownames(full_list_efit_IvsK), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano I vs K (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_IvsM,
                lab = rownames(full_list_efit_IvsM), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano I vs M (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

EnhancedVolcano(full_list_efit_MvsK,
                lab = rownames(full_list_efit_MvsK), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-4, 4), 
                ylim = c(0, 7),
                FCcutoff = 2,  
                title ="Vulcano M vs K (HCA)",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)


#Heatmap
  p_above_FC <- rownames(full_list_efit_CvsM)[(full_list_efit_CvsM$adj.P.Val <= 0.05)&((full_list_efit_CvsM$logFC >= 2)|(full_list_efit_CvsM$logFC <= -2))]
  length(p_above_FC) 
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_CvsM)[(full_list_efit_CvsM$adj.P.Val <= 0.05)&((full_list_efit_CvsM$logFC >= 2)|(full_list_efit_CvsM$logFC <= -2))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap C vs M")

  p_above_FC <- rownames(full_list_efit_IvsC)[(full_list_efit_IvsC$adj.P.Val <= 0.05)&((full_list_efit_IvsC$logFC >= 1)|(full_list_efit_IvsC$logFC <= -1))]
  length(p_above_FC)
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_IvsC)[(full_list_efit_IvsC$adj.P.Val <= 0.05)&((full_list_efit_IvsC$logFC >= 1)|(full_list_efit_IvsC$logFC <= -1))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap I vs C LogFC(1;-1)")

  p_above_FC <- rownames(full_list_efit_CvsK)[(full_list_efit_CvsK$adj.P.Val <= 0.05)&((full_list_efit_CvsK$logFC >= 1)|(full_list_efit_CvsK$logFC <= -1))]
  length(p_above_FC) 
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_CvsK)[(full_list_efit_CvsK$adj.P.Val <= 0.05)&((full_list_efit_CvsK$logFC >= 2)|(full_list_efit_CvsK$logFC <= -2))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap C vs K")

  p_above_FC <- rownames(full_list_efit_IvsK)[(full_list_efit_IvsK$adj.P.Val <= 0.05)&((full_list_efit_IvsK$logFC >= 1)|(full_list_efit_IvsK$logFC <= -1))]
  length(p_above_FC) 
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_IvsK)[(full_list_efit_IvsK$adj.P.Val <= 0.05)&((full_list_efit_IvsK$logFC >= 2)|(full_list_efit_IvsK$logFC <= -2))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap I vs K")

  p_above_FC <- rownames(full_list_efit_IvsM)[(full_list_efit_IvsM$adj.P.Val <= 0.05)&((full_list_efit_IvsM$logFC >= 1)|(full_list_efit_IvsM$logFC <= -1))]
  length(p_above_FC) 
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_IvsM)[(full_list_efit_IvsM$adj.P.Val <= 0.05)&((full_list_efit_IvsM$logFC >= 2)|(full_list_efit_IvsM$logFC <= -2))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap I vs M")

  p_above_FC <- rownames(full_list_efit_MvsK)[(full_list_efit_MvsK$adj.P.Val <= 0.05)&((full_list_efit_MvsK$logFC >= 1)|(full_list_efit_MvsK$logFC <= -1))]
  length(p_above_FC)
  p_above_FC_top_for_heatmap <- rownames(full_list_efit_MvsK)[(full_list_efit_MvsK$adj.P.Val <= 0.05)&((full_list_efit_MvsK$logFC >= 2)|(full_list_efit_MvsK$logFC <= -2))]
  length(p_above_FC_top_for_heatmap)
pheatmap(dat_norm[p_above_FC_top_for_heatmap,], annotation_col = fact[,c(2,1)], main = "Heatmap M vs K")

