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

fact <- data.frame(read_excel("sample_info_131-153.xlsx")) 
  
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
  
  fact$Original.index_2 <- paste(fact$Differentiation, fact$Original.index, sep = "_") #workname in registry
  fact$Original.index_2 <- as.factor(fact$Original.index_2)
  
  rownames(fact) <- fact$Original.index_2
  str(fact)
  
dat <- data.frame(read.delim("proteins_131-153.tsv")) 
  dat[dat==0] <- NA 

  dat1 <- dat[,c(1,103:124)]
  name_gene <- dat[,c(2,4)]
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
      # HCA vs HIT
    {  
    df <- data.frame(part=0:1000)
    c <- matrix(,ncol = 1001, nrow=2)
    for (i in 0:1000) 
      {
      HIT <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HIT"))])) >= i/1000), ]
      HCA <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HCA"))])) >= i/1000), ]
      c[1,i+1] <- nrow(HIT)
      c[2,i+1] <- nrow(HCA)
      }
    df['HIT'] <- c[1,]
    df['HCA'] <- c[2,]
    ggplot(data = df) +
      geom_smooth(mapping = aes(x = part, y = HIT, colour = "HIT")) +             
      geom_line(mapping = aes(x = part, y = HIT, colour = "HIT")) +    
      geom_smooth(mapping = aes(x = part, y = HCA, colour = "HCA")) + 
      geom_line(mapping = aes(x = part, y = HCA, colour = "HCA")) +
      geom_vline(xintercept = 600) +
      xlab("Decline NA (‰)") + ylab("Detected proteins") +
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      scale_y_continuous(breaks = seq(0, 3000, 500)) +
      labs(subtitle="Dependence of the number of detectable proteins on decline NA (all HCA + all HIT)")
    }

      # all C vs K vs M vs I
    {  
    df <- data.frame(part=0:1000)
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
      geom_smooth(mapping = aes(x = part, y = I, colour = "I"), se = F) + 
      geom_line(mapping = aes(x = part, y = I, colour = "I")) +
      geom_smooth(mapping = aes(x = part, y = C, colour = "C"), se = F) + 
      geom_line(mapping = aes(x = part, y = C, colour = "C")) +
      geom_smooth(mapping = aes(x = part, y = M, colour = "M"), se = F) + 
      geom_line(mapping = aes(x = part, y = M, colour = "M")) +
      geom_smooth(mapping = aes(x = part, y = K, colour = "K"), se = F) + 
      geom_line(mapping = aes(x = part, y = K, colour = "K")) +
      geom_vline(xintercept = 600) +
      xlab("Decline NA (‰)") + ylab("Detected proteins") +
      scale_x_continuous(breaks = seq(0, 1000, 100)) +
      scale_y_continuous(breaks = seq(0, 3000, 500)) +
      labs(subtitle="Dependence of the number of detectable proteins on decline NA (all HCA + all HIT) by groups")
    }

      # HIT+HCA: C vs K vs M vs I
    {  
      df <- data.frame(part=0:1000)
      c <- matrix(,ncol = 1001, nrow=8)
      for (i in 0:1000) 
      {
        I <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_I"))])) >= i/1000), ]
        C <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_C"))])) >= i/1000), ]
        M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_M"))])) >= i/1000), ]
        K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HCA_K"))])) >= i/1000), ]
        II <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_I"))])) >= i/1000), ]
        CC <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_C"))])) >= i/1000), ]
        MM <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_M"))])) >= i/1000), ]
        KK <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group_2=="HIT_K"))])) >= i/1000), ]
        c[1,i+1] <- nrow(I)
        c[2,i+1] <- nrow(C)
        c[3,i+1] <- nrow(M)
        c[4,i+1] <- nrow(K)
        c[5,i+1] <- nrow(II)
        c[6,i+1] <- nrow(CC)
        c[7,i+1] <- nrow(MM)
        c[8,i+1] <- nrow(KK)
      }
      df['HCA_I'] <- c[1,]
      df['HCA_C'] <- c[2,]
      df['HCA_M'] <- c[3,]
      df['HCA_K'] <- c[4,]
      df['HIT_I'] <- c[5,]
      df['HIT_C'] <- c[6,]
      df['HIT_M'] <- c[7,]
      df['HIT_K'] <- c[8,]
      ggplot(data = df) +
        #geom_smooth(mapping = aes(x = part, y = HCA_I, colour = "HCA_I"), se = F) + 
        geom_line(mapping = aes(x = part, y = HCA_I, colour = "HCA_I")) +
        #geom_smooth(mapping = aes(x = part, y = HCA_C, colour = "HCA_C"), se = F) + 
        geom_line(mapping = aes(x = part, y = HCA_C, colour = "HCA_C")) +
        #geom_smooth(mapping = aes(x = part, y = HCA_M, colour = "HCA_M"), se = F) + 
        geom_line(mapping = aes(x = part, y = HCA_M, colour = "HCA_M")) +
        #geom_smooth(mapping = aes(x = part, y = HCA_K, colour = "HCA_K"), se = F) + 
        geom_line(mapping = aes(x = part, y = HCA_K, colour = "HCA_K")) +
        #geom_smooth(mapping = aes(x = part, y = HIT_I, colour = "HIT_I"), se = F) + 
        geom_line(mapping = aes(x = part, y = HIT_I, colour = "HIT_I")) +
        #geom_smooth(mapping = aes(x = part, y = HIT_C, colour = "HIT_C"), se = F) + 
        geom_line(mapping = aes(x = part, y = HIT_C, colour = "HIT_C")) +
        #geom_smooth(mapping = aes(x = part, y = HIT_M, colour = "HIT_M"), se = F) + 
        geom_line(mapping = aes(x = part, y = HIT_M, colour = "HIT_M")) +
        #geom_smooth(mapping = aes(x = part, y = HIT_K, colour = "HIT_K"), se = F) + 
        geom_line(mapping = aes(x = part, y = HIT_K, colour = "HIT_K")) +
        geom_vline(xintercept = 600) +
        xlab("Decline NA (‰)") + ylab("Detected proteins") +
        scale_x_continuous(breaks = seq(0, 1000, 100)) +
        scale_y_continuous(breaks = seq(0, 3000, 500)) +
        labs(subtitle="Dependence of the number of detectable proteins on decline NA (all HCA + all HIT)")
    }

  #filter out proteins with too much NA
  HIT <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HIT"))])) >= 0.71), ]
  HCA <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,Differentiation=="HCA"))])) >= 0.71), ]
  C <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="C"))])) >= 0.5), ]
  K <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="K"))])) >= 0.5), ]
  M <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="M"))])) >= 0.5), ]
  I <- dat1[which(rowMeans(!is.na(dat1[,rownames(subset(fact,group=="I"))])) >= 0.5), ]
  
  vennn <- list(HCA = rownames(HCA), HIT = rownames(HIT))
  vennn2 <- list(C = rownames(C), K = rownames(K), M = rownames(M), I = rownames(I))
  vennn3 <- list(C = rownames(C), K = rownames(K), I = rownames(I))

ggvenn(vennn, 
       fill_color = c("#0073C2FF", "#CD534CFF"),
       stroke_size = 0.5, set_name_size = 8, text_size = 7,)

ggvenn(vennn2, 
       fill_color = c("#0073C2FF", "#CD534CFF", "#009E73", "#F0E442"),
       stroke_size = 0.5, set_name_size = 8, text_size = 5,)

ggvenn(vennn3, 
       fill_color = c("#0073C2FF", "#CD534CFF", "#009E73"),
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
  
  HCA_spec <- Setdiff(vennn[c("HCA")], vennn[c("HIT")])
  HIT_spec <- Setdiff(vennn[c("HIT")], vennn[c("HCA")])
  C_spec <- Setdiff(vennn2[c("C")], vennn2[c("K", 'M', 'I')])
  M_spec <- Setdiff(vennn2[c("M")], vennn2[c("K", 'C', 'I')])
  I_spec <- Setdiff(vennn2[c("I")], vennn2[c("K", 'M', 'C')])
  K_spec <- Setdiff(vennn2[c("K")], vennn2[c("C", 'M', 'I')])
  
  # write.table(HCA_spec, "HCA_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(HIT_spec, "HIT_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(C_spec, "C_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(I_spec, "I_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(K_spec, "K_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  # write.table(M_spec, "M_spec_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

#Quantitative analysis

  dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.7),] #filtration data according to the above algorithm
  mean(complete.cases(dat2))
  NAsums <- data.frame(colSums(is.na(dat2))) #detection NA in experimental sample
  NAsums 
  # write.table(NAsums, "NA_sums_131-153.txt", append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  str(dat2)

#Imputation
  tdat <- t(dat2)
  dat_knn1 <- impute.knn(tdat, k = 5)
  dat_knn <- t(dat_knn1$data)
  mean(complete.cases(dat_knn))

#Structure of row data
  pal <- brewer.pal(n = 9, name = "Set1")
  
  cols <- pal[fact$Differentiation]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data all HCA+HIT")
  legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
  
  cols <- pal[fact$group_2]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data all HCA+HIT")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)
  
  colSums(dat_knn)
  length(dat_knn[dat_knn == 0]) 

#Logarithm of data
  dat_log <- log2(dat_knn+1) 
  head(dat_log)
  mean(complete.cases(dat_log))
  
  cols <- pal[fact$Differentiation]
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data all HCA+HIT")
  legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
  
  cols <- pal[fact$group_2]
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data all HCA+HIT")
  legend("topright", levels(fact$group_2), fill = pal, bty = "n", xpd = T)

#Normalization of data
  dat_norm <- normalizeQuantiles(dat_log) 
  head(dat_norm)

  cols <- pal[fact$Differentiation]
boxplot(dat_norm, col = cols, main = "Normalized data all HCA+HIT")
  legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
  
  cols <- pal[fact$group_2]
boxplot(dat_norm, col = cols, main = "Normalized data all HCA+HIT")
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

maplot(dat_log[, rownames(fact)[fact$Differentiation == "HIT"]], dat_log[, rownames(fact)[fact$Differentiation == "HCA"]], main = "Log-expression data all HCA+HIT")
maplot(dat_norm[, rownames(fact)[fact$Differentiation == "HIT"]], dat_norm[, rownames(fact)[fact$Differentiation == "HCA"]], main = "Normalized data all HCA+HIT")

#MeanSd
meanSdPlot(as.matrix(dat_log))
meanSdPlot(as.matrix(dat_norm))

#Heatmap
aheatmap(cor(dat_norm), 
         color = c("#00bfff", '#240935',"#e600e6"), 
         annCol = fact$Differentiation, 
         fontsize = 10)

#Principle components analysis
  dat_pca <- pca(t(dat_norm), ncomp = 10, center = TRUE)
  dat_pca

plot(dat_pca)

plotIndiv(dat_pca, comp = c(1, 2),
          ind.names = FALSE, 
          group = fact$group,
          legend = TRUE,
          ellipse = TRUE,
          title = 'PCA HCA vs HIT by group') 

plotIndiv(dat_pca, comp = c(1, 2),
          ind.names = FALSE, 
          group = fact$group_2, 
          legend = TRUE,
          ellipse = TRUE,
          title = 'PCA HCA vs HIT by undergroup',
          ellipse.level = 0.75)

plotIndiv(dat_pca, comp = c(1, 2), 
          ind.names = FALSE, 
          group = fact$Differentiation, 
          legend = TRUE, 
          ellipse = TRUE,
          title = 'PCA HCA vs HIT', 
          ellipse.level = 0.95) 

#Limma - differentially expressed proteins
  X <- model.matrix(~ fact$Differentiation)
  X
  
  fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)
  
  efit <- eBayes(fit)
  
  topTable(efit, coef = 2)
  numGenes <- length(dat_norm)
  full_list_efit <- topTable(efit, number = length(dat_norm[,1]))
  # write.csv(full_list_efit,'Dif_expr_HCA_vs_HIT.csv')
  head(full_list_efit)
  
  full_list_efit$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", rownames(full_list_efit)), column = "SYMBOL", keytype = "UNIPROT")
  sum(is.na(full_list_efit$Gene.name)) 
  full_list_efit$Gene.name[is.na(full_list_efit$Gene.name)] <- sub("\\|.*$", "",rownames(full_list_efit)[is.na(full_list_efit$Gene.name)])
  full_list_efit['Name.of.protein'] <- sub('\\|.*$', '', rownames(full_list_efit))
  rownames(full_list_efit) <- full_list_efit$Gene.name
  
  # write.xlsx(full_list_efit, "HCA_vs_HIT_diff_proteins.xlsx",append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
  
#Vulcano plot
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit), 
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-3, 3), 
                ylim = c(0, 8),
                FCcutoff = 1,  
                title ="Vulcano all HCA+HIT",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)

#Heatmap

  #rename values for building modified heatmap legend
  full_list_efit2 <- full_list_efit[order(full_list_efit$Name.of.protein),]
  clusters_col <- c('HCA_K2', 'HCA_K4', 'HCA_M2', 'HCA_M6', 'HIT_K2', 'HIT_K4', 'HIT_K6', 'HIT_M2', 'HIT_M4', 'HIT_M6', 'HCA_C2', 'HCA_C4', 'HCA_C6', 'HCA_I2', 'HCA_I4', 'HCA_I6', 'HIT_C2', 'HIT_C4', 'HIT_C6', 'HIT_I2', 'HIT_I4', 'HIT_I6')
  print(colnames(dat2))
  print(colnames(dat_norm))
  dat2_2 <- dat2[,clusters_col]
  dat_norm2 <- dat_norm[,clusters_col]
  print(colnames(dat2_2))
  print(colnames(dat_norm2))
  full_list_efit3 <- full_list_efit
  rownames(full_list_efit3) <- paste(full_list_efit$Name.of.protein, full_list_efit$Gene.name, sep = '--')
  rownames(dat2_2) <- paste(full_list_efit2$Name.of.protein,full_list_efit2$Gene.name, sep = '--')
  rownames(dat_norm2) <- paste(full_list_efit2$Name.of.protein,full_list_efit2$Gene.name, sep = '--')
  
  #selected proteins for heatmap visualization. They be selected on the basis of average expression 
  prot_for_heatmap <- rownames(full_list_efit3[c('P13987','P14174','Q03135','P08962','Q14254',          #top5 of overexpression
                                                 'P09486','P98160','P02751','P01033','P07942','Q16363', #basement membrane
                                                 'P11047','Q9Y4K0','P14543','P35555','O00468','P08572',
                                                 'P12109','Q92626',
                                                 'P07602','P07996','P29279','Q12805','P04275','Q16270', #extracellular matrix
                                                 'P07900','P02765','P16035','P08253',
                                                 'P27797','O00622','P10909','Q13201','P21810','P55001','P51858','P14625','Q9GZM7',
                                                 'P50454','P05997','Q14767','P19823','Q15063',
                                                 'P07711','O60568','Q7Z7G0'),])
  length(prot_for_heatmap)

pheatmap(dat_norm2[prot_for_heatmap,], #main heatmap
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm2[prot_for_heatmap,], #main heatmap in green-red gradient (experimentally)
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00ff00",'black','#ff0000'))(100))


pheatmap(dat2_2[prot_for_heatmap,], #main heatmap which was builded on base area of peak in spectral analyse (sqr(experimentally))
         annotation_col = fact[,c(2,5)],
         cutree_cols = 1,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         cellwidth = 12,
         cellheight = 10,
         breaks = seq(0, 2000000, length.out = 100),
         border_color = "black",
         color = colorRampPalette(c('#005aeb', '#240935', '#800080', '#b300b3',"#e600e6"))(100))


  pph <- rownames(head(full_list_efit3[order(full_list_efit3$logFC,decreasing = T),],n=50)) #rownames of top50 overexpression proteins
  length(pph)
  
pheatmap(dat_norm2[pph,], #heatmap of top50 overexpression proteins (experimentally)
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))


  #the spread of legend values
  min(dat_norm2[prot_for_heatmap,])
  max(dat_norm2[prot_for_heatmap,])

#Unique proteins and preparing table
  
  #load proteins which were analysed in the last publication (PRIDE: PXD038017)
  dat_prot_HCA <- data.frame(read.csv('proteome_HCA.csv'))
    dat_prot_HCA$Protein <- sub("\\|.*$", "",dat_prot_HCA$Accession)
    dat_prot_HCA$UNIPROT <- sub(".*\\|", "",dat_prot_HCA$Accession)
    dat_prot_HCA$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", dat_prot_HCA$Accession), column = "SYMBOL", keytype = "UNIPROT")
  dat_prot_HIT <- data.frame(read.csv('proteome_HIT.csv'))
    dat_prot_HIT$Protein <- sub("\\|.*$", "",dat_prot_HIT$Accession)
    dat_prot_HIT$UNIPROT <- sub(".*\\|", "",dat_prot_HIT$Accession)
    dat_prot_HIT$Gene.name <- mapIds(org.Hs.eg.db, keys = sub("\\|.*$", "", dat_prot_HIT$Accession), column = "SYMBOL", keytype = "UNIPROT")
  
  #Areas and Unique proteins
  
  Uniq_prot_secretome <- dat[,c(4,2,3,11,103:124,15:36)] # All secretome proteins for HCA+HIT from MSfragger-output
  
  Uniq_prot_proteome <- dat_prot_HCA[,c(49,47,48,5)] # All proteome proteins for HCA+HIT from PXD038017
    Uniq_prot_proteome_uniq_pep_HIT <- dat_prot_HIT[,c(49,47,48,5)] #uniq peptides of HITs
    Uniq_prot_proteome_areas_HCA <- dat_prot_HCA[,c(49,47,48,24:29,18:23)] #areas of HCAs
    Uniq_prot_proteome_areas_HIT <- dat_prot_HIT[,c(49,47,48,24:29,18:23)] #areas of HITs
    Uniq_prot_proteome_count_HCA <- dat_prot_HCA[,c(49,47,48,12:17,6:11)] #counts of HCAs
    Uniq_prot_proteome_count_HIT <- dat_prot_HIT[,c(49,47,48,12:17,6:11)] #counts of HITs
    Uniq_prot_proteome <- merge(Uniq_prot_proteome, Uniq_prot_proteome_uniq_pep_HIT, by = c('Gene.name','Protein','UNIPROT'), all = T)
    Uniq_prot_proteome <- merge(Uniq_prot_proteome, Uniq_prot_proteome_areas_HCA, by = c('Gene.name','Protein','UNIPROT'), all = T)
    Uniq_prot_proteome <- merge(Uniq_prot_proteome, Uniq_prot_proteome_areas_HIT, by = c('Gene.name','Protein','UNIPROT'), all = T)
    Uniq_prot_proteome <- merge(Uniq_prot_proteome, Uniq_prot_proteome_count_HCA, by = c('Gene.name','Protein','UNIPROT'), all = T)
    Uniq_prot_proteome <- merge(Uniq_prot_proteome, Uniq_prot_proteome_count_HIT, by = c('Gene.name','Protein','UNIPROT'), all = T)
  
  Uniq_prot_secretome <- Uniq_prot_secretome[which((rowMeans(!is.na(Uniq_prot_secretome[,5:26])) >0.71)),] # Secretome proteins for HCA+HIT which intensity 1 more
  Uniq_prot_secretome <- Uniq_prot_secretome[which(Uniq_prot_secretome$Combined.Total.Peptides > 1),] #delete proteins who have uniq peptides <2
  Uniq_prot_proteome <- Uniq_prot_proteome[which((rowMeans(!is.na(Uniq_prot_proteome[,6:29])) >0.9)),] # Proteome proteins for HCA+HIT which intensity 1 more
  Uniq_prot_proteome <- Uniq_prot_proteome[which((Uniq_prot_proteome$Coverage.....x > 1) | (Uniq_prot_proteome$Coverage.....y > 1)),] #delete proteins who have uniq peptides <2 in HCA+HIT
  
  colnames(Uniq_prot_secretome) <- c('Gene name',
                                     'Protein',
                                     'UNIPROT',
                                     'Peptides',
                                     'Area HCAEC CPP-P 2','Area HCAEC CPP-P 4','Area HCAEC CPP-P 6','Area HCAEC CPP-S 2','Area HCAEC CPP-S 4','Area HCAEC CPP-S 6','Area HCAEC DPBS 2','Area HCAEC DPBS 4','Area HCAEC MPP 2','Area HCAEC MPP 6',
                                     'Area HITAEC CPP-P 2','Area HITAEC CPP-P 4','Area HITAEC CPP-P 6','Area HITAEC CPP-S 2','Area HITAEC CPP-S 4','Area HITAEC CPP-S 6','Area HITAEC DPBS 2','Area HITAEC DPBS 4','Area HITAEC DPBS 6','Area HITAEC MPP 2','Area HITAEC MPP 4','Area HITAEC MPP 6',
                                     'Count HCAEC CPP-P 2','Count HCAEC CPP-P 4','Count HCAEC CPP-P 6','Count HCAEC CPP-S 2','Count HCAEC CPP-S 4','Count HCAEC CPP-S 6','Count HCAEC DPBS 2','Count HCAEC DPBS 4','Count HCAEC MPP 2','Count HCAEC MPP 6',
                                     'Count HITAEC CPP-P 2','Count HITAEC CPP-P 4','Count HITAEC CPP-P 6','Count HITAEC CPP-S 2','Count HITAEC CPP-S 4','Count HITAEC CPP-S 6','Count HITAEC DPBS 2','Count HITAEC DPBS 4','Count HITAEC DPBS 6','Count HITAEC MPP 2','Count HITAEC MPP 4','Count HITAEC MPP 6')
  colnames(Uniq_prot_proteome) <- c('Gene name',
                                    'Protein',
                                    'UNIPROT',
                                    'Peptides HCAEC','Peptides HITAEC',
                                    'Area HCAEC CPP-P 2','Area HCAEC CPP-P 4','Area HCAEC CPP-P 6','Area HCAEC CPP-S 2','Area HCAEC CPP-S 4','Area HCAEC CPP-S 6','Area HCAEC DPBS 2','Area HCAEC DPBS 4','Area HCAEC DPBS 6','Area HCAEC MPP 2','Area HCAEC MPP 4','Area HCAEC MPP 6',
                                    'Area HITAEC CPP-P 2','Area HITAEC CPP-P 4','Area HITAEC CPP-P 6','Area HITAEC CPP-S 2','Area HITAEC CPP-S 4','Area HITAEC CPP-S 6','Area HITAEC DPBS 2','Area HITAEC DPBS 4','Area HITAEC DPBS 6','Area HITAEC MPP 2','Area HITAEC MPP 4','Area HITAEC MPP 6',
                                    'Count HCAEC CPP-P 2','Count HCAEC CPP-P 4','Count HCAEC CPP-P 6','Count HCAEC CPP-S 2','Count HCAEC CPP-S 4','Count HCAEC CPP-S 6','Count HCAEC DPBS 2','Count HCAEC DPBS 4','Count HCAEC DPBS 6','Count HCAEC MPP 2','Count HCAEC MPP 4','Count HCAEC MPP 6',
                                    'Count HITAEC CPP-P 2','Count HITAEC CPP-P 4','Count HITAEC CPP-P 6','Count HITAEC CPP-S 2','Count HITAEC CPP-S 4','Count HITAEC CPP-S 6','Count HITAEC DPBS 2','Count HITAEC DPBS 4','Count HITAEC DPBS 6','Count HITAEC MPP 2','Count HITAEC MPP 4','Count HITAEC MPP 6')
  
  Uniq_prot_secretome_exp_HIT <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,15:20])) == 6) & #HIT exp
                                   (rowSums(!is.na(Uniq_prot_secretome[,21:26])) == 0), -c(5:14,27:36)] #HITAEC Ctrl and delete HCA
  
  Uniq_prot_secretome_ctrl_HIT <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,21:26])) == 6) & 
                                    (rowSums(!is.na(Uniq_prot_secretome[,15:20])) == 0), -c(4:13,27:36)]
  
  Uniq_prot_secretome_exp_HCA <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,5:10])) == 6) & 
                                   (rowSums(!is.na(Uniq_prot_secretome[,11:14])) == 0), -c(15:26,37:48)] 
  
  Uniq_prot_secretome_ctrl_HCA <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,11:14])) == 4) & 
                                    (rowSums(!is.na(Uniq_prot_secretome[,5:10])) == 0), -c(15:26,37:48)]
  
  Uniq_prot_secretome_HCA_ctrl_vs_HIT_ctrl <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,11:14])) == 4) & 
                                   (rowSums(!is.na(Uniq_prot_secretome[,21:26])) == 0), -c(5:10,15:20,27:32,37:42)] 
  
  Uniq_prot_secretome_HIT_ctrl_vs_HCA_ctrl <- Uniq_prot_secretome[(rowSums(!is.na(Uniq_prot_secretome[,21:25])) == 6) & 
                                    (rowSums(!is.na(Uniq_prot_secretome[,11:14])) == 0), -c(5:10,15:20,27:32,37:42)]
  
  #preparing data for export to excel-book and correlation analysis in Prism
  
  Uniq_prot_secretome[is.na(Uniq_prot_secretome)] <- 0 # NAs to 0 for correct correlation analyse
  Uniq_prot_proteome[is.na(Uniq_prot_proteome)] <- 0
  
  {
  Uniq_prot_areas_proteome <- Uniq_prot_proteome[1:3]
    Uniq_prot_areas_proteome$'MeanArea HCAEC CPP-P' <- rowMeans(Uniq_prot_proteome[,c(6:8)])
    Uniq_prot_areas_proteome$'MeanArea HCAEC CPP-S' <- rowMeans(Uniq_prot_proteome[,c(9:11)])
    Uniq_prot_areas_proteome[['MeanArea HCAEC (CPP-P + CPP-S)']] <- rowMeans(Uniq_prot_proteome[,c(6:11)])
    Uniq_prot_areas_proteome$'MeanArea HCAEC DPBS' <- rowMeans(Uniq_prot_proteome[,c(12:14)])
    Uniq_prot_areas_proteome$'MeanArea HCAEC MPP' <- rowMeans(Uniq_prot_proteome[,c(15:17)])
    Uniq_prot_areas_proteome[['MeanArea HCAEC (DPBS + MPP)']] <- rowMeans(Uniq_prot_proteome[,c(12:17)])
    Uniq_prot_areas_proteome$'MeanArea HITAEC CPP-P' <- rowMeans(Uniq_prot_proteome[,c(18:20)])
    Uniq_prot_areas_proteome$'MeanArea HITAEC CPP-S' <- rowMeans(Uniq_prot_proteome[,c(21:23)])
    Uniq_prot_areas_proteome[['MeanArea HITAEC (CPP-P + CPP-S)']] <- rowMeans(Uniq_prot_proteome[,c(18:23)])
    Uniq_prot_areas_proteome$'MeanArea HITAEC DPBS' <- rowMeans(Uniq_prot_proteome[,c(24:26)])
    Uniq_prot_areas_proteome$'MeanArea HITAEC MPP' <- rowMeans(Uniq_prot_proteome[,c(27:29)])
    Uniq_prot_areas_proteome[['MeanArea HITAEC (DPBS + MPP)']] <- rowMeans(Uniq_prot_proteome[,c(24:29)])
  }
  
  {
  Uniq_prot_areas_secretome <- Uniq_prot_secretome[1:3]
    Uniq_prot_areas_secretome$'MeanArea HCAEC CPP-P' <- rowMeans(Uniq_prot_secretome[,c(5:7)])
    Uniq_prot_areas_secretome$'MeanArea HCAEC CPP-S' <- rowMeans(Uniq_prot_secretome[,c(8:10)])
    Uniq_prot_areas_secretome[['MeanArea HCAEC (CPP-P + CPP-S)']] <- rowMeans(Uniq_prot_secretome[,c(5:10)])
    Uniq_prot_areas_secretome$'MeanArea HCAEC DPBS' <- rowMeans(Uniq_prot_secretome[,c(11:12)])
    Uniq_prot_areas_secretome$'MeanArea HCAEC MPP' <- rowMeans(Uniq_prot_secretome[,c(13:14)])
    Uniq_prot_areas_secretome[['MeanArea HCAEC (DPBS + MPP)']] <- rowMeans(Uniq_prot_secretome[,c(11:14)])
    Uniq_prot_areas_secretome$'MeanArea HITAEC CPP-P' <- rowMeans(Uniq_prot_secretome[,c(15:17)])
    Uniq_prot_areas_secretome$'MeanArea HITAEC CPP-S' <- rowMeans(Uniq_prot_secretome[,c(18:20)])
    Uniq_prot_areas_secretome[['MeanArea HITAEC (CPP-P + CPP-S)']] <- rowMeans(Uniq_prot_secretome[,c(15:20)])
    Uniq_prot_areas_secretome$'MeanArea HITAEC DPBS' <- rowMeans(Uniq_prot_secretome[,c(21:23)])
    Uniq_prot_areas_secretome$'MeanArea HITAEC MPP' <- rowMeans(Uniq_prot_secretome[,c(24:26)])
    Uniq_prot_areas_secretome[['MeanArea HITAEC (DPBS + MPP)']] <-rowMeans(Uniq_prot_secretome[,c(21:26)])
  }
    
  prot_heatmap <- c('P13987','P14174','Q03135','P08962','Q14254',          #over
                    'P09486','P98160','P02751','P01033','P07942','Q16363', #membrane
                    'P11047','Q9Y4K0','P14543','P35555','O00468','P08572',
                    'P12109','Q92626',
                    'P07602','P07996','P29279','Q12805','P04275','Q16270', #matrix
                    'P07900','P02765','P16035','P08253',
                    'P27797','O00622','P10909','Q13201','P21810','P55001','P51858','P14625','Q9GZM7',
                    'P50454','P05997','Q14767','P19823','Q15063',
                    'P07711','O60568','Q7Z7G0')
  length(prot_heatmap)
    
    ctrl_heatmap_prot <- Uniq_prot_proteome[,c(1:3,6:29)][Uniq_prot_proteome$Protein %in% prot_heatmap,]
    ctrl_heatmap_secr <- Uniq_prot_secretome[,c(1:3,5:26)][Uniq_prot_secretome$Protein %in% prot_heatmap,]
  all_heatmap <- merge(ctrl_heatmap_secr,ctrl_heatmap_prot, by = c('Gene name','Protein','UNIPROT'), all.x = T)
  all_heatmap[is.na(all_heatmap)] <- 0
  
  all_protein <- merge(Uniq_prot_secretome[,c(1:3,5:26)], Uniq_prot_proteome[,c(1:3,6:29)], by = c('Gene name','Protein','UNIPROT'))
  
  {
  all_protein_means <- all_protein[,c(1:3)]
    all_protein_means$'Secr meanArea HCAEC CPP-P' <- rowMeans(all_protein[,c(4:6)])
    all_protein_means$'Secr meanArea HCAEC CPP-S' <- rowMeans(all_protein[,c(7:9)])
    all_protein_means$'Secr meanArea HCAEC (CPP-P + CPP-S)' <- rowMeans(all_protein[,c(4:9)])
    all_protein_means$'Secr meanArea HCAEC DPBS' <- rowMeans(all_protein[,c(10:11)])
    all_protein_means$'Secr meanArea HCAEC MPP' <- rowMeans(all_protein[,c(12:13)])
    all_protein_means$'Secr meanArea HCAEC (DPBS + MPP)' <- rowMeans(all_protein[,c(10:13)])
    all_protein_means$'Secr meanArea HITAEC CPP-P' <- rowMeans(all_protein[,c(14:16)])
    all_protein_means$'Secr meanArea HITAEC CPP-S' <- rowMeans(all_protein[,c(17:19)])
    all_protein_means$'Secr meanArea HITAEC (CPP-P + CPP-S)' <- rowMeans(all_protein[,c(14:19)])
    all_protein_means$'Secr meanArea HITAEC DPBS' <- rowMeans(all_protein[,c(20:22)])
    all_protein_means$'Secr meanArea HITAEC MPP' <- rowMeans(all_protein[,c(23:25)])
    all_protein_means$'Secr meanArea HITAEC (DPBS + MPP)' <- rowMeans(all_protein[,c(20:25)])
    all_protein_means$'Secr meanArea' <- rowMeans(all_protein[,c(4:25)])
    all_protein_means$'Prot meanArea HCAEC CPP-P' <- rowMeans(all_protein[,c(26:28)])
    all_protein_means$'Prot meanArea HCAEC CPP-S' <- rowMeans(all_protein[,c(29:31)])
    all_protein_means$'Prot meanArea HCAEC (CPP-P + CPP-S)' <- rowMeans(all_protein[,c(26:31)])
    all_protein_means$'Prot meanArea HCAEC DPBS' <- rowMeans(all_protein[,c(32:34)])
    all_protein_means$'Prot meanArea HCAEC MPP' <- rowMeans(all_protein[,c(35:37)])
    all_protein_means$'Prot meanArea HCAEC (DPBS + MPP)' <- rowMeans(all_protein[,c(32:37)])
    all_protein_means$'Prot meanArea HITAEC CPP-P' <- rowMeans(all_protein[,c(38:40)])
    all_protein_means$'Prot meanArea HITAEC CPP-S' <- rowMeans(all_protein[,c(41:43)])
    all_protein_means$'Prot meanArea HITAEC (CPP-P + CPP-S)' <- rowMeans(all_protein[,c(38:43)])
    all_protein_means$'Prot meanArea HITAEC DPBS' <- rowMeans(all_protein[,c(44:46)])
    all_protein_means$'Prot meanArea HITAEC MPP' <- rowMeans(all_protein[,c(47:49)])
    all_protein_means$'Prot meanArea HITAEC (DPBS + MPP)' <- rowMeans(all_protein[,c(44:49)])
    all_protein_means$'Prot meanArea' <- rowMeans(all_protein[,c(26:49)])
  }
    
  {  
  all_heatmap_means <- all_heatmap[,c(1:3)]
    all_heatmap_means$'Secr meanArea HCAEC (CPP-P + CPP-S)' <- rowMeans(all_heatmap[,c(4:9)])
    all_heatmap_means$'Secr meanArea HCAEC (DPBS + MPP)' <- rowMeans(all_heatmap[,c(10:13)])
    all_heatmap_means$'Secr meanArea HITAEC (CPP-P + CPP-S)' <- rowMeans(all_heatmap[,c(14:19)])
    all_heatmap_means$'Secr meanArea HITAEC (DPBS + MPP)' <- rowMeans(all_heatmap[,c(20:25)])
    all_heatmap_means$'Secr meanArea' <- rowMeans(all_heatmap[,c(4:25)])
    all_heatmap_means$'Prot meanArea HCAEC (CPP-P + CPP-S)' <- rowMeans(all_heatmap[,c(26:31)])
    all_heatmap_means$'Prot meanArea HCAEC (DPBS + MPP)' <- rowMeans(all_heatmap[,c(32:37)])
    all_heatmap_means$'Prot meanArea HITAEC (CPP-P + CPP-S)' <- rowMeans(all_heatmap[,c(38:43)])
    all_heatmap_means$'Prot meanArea HITAEC (DPBS + MPP)' <- rowMeans(all_heatmap[,c(44:49)])
    all_heatmap_means$'Prot meanArea' <- rowMeans(all_heatmap[,c(26:49)])
  }
  
  
  #Average expression
     
  all_secr_AE <- data.frame(dat_norm2)
    all_secr_AE$'Gene name' <- sub(".*\\-", "", rownames(all_secr_AE))
    all_secr_AE$'Protein' <- sub("\\-.*$", "",rownames(all_secr_AE))
    all_secr_AE <- all_secr_AE[,c(23,24,11:16,1:4,17:22,5:10)]
    
    colnames(all_secr_AE) <- c('Gene name',
                                       'Protein',
                                       'AE HCAEC CPP-P 2','AE HCAEC CPP-P 4','AE HCAEC CPP-P 6','AE HCAEC CPP-S 2','AE HCAEC CPP-S 4','AE HCAEC CPP-S 6','AE HCAEC DPBS 2','AE HCAEC DPBS 4','AE HCAEC MPP 2','AE HCAEC MPP 6',
                                       'AE HITAEC CPP-P 2','AE HITAEC CPP-P 4','AE HITAEC CPP-P 6','AE HITAEC CPP-S 2','AE HITAEC CPP-S 4','AE HITAEC CPP-S 6','AE HITAEC DPBS 2','AE HITAEC DPBS 4','AE HITAEC DPBS 6','AE HITAEC MPP 2','AE HITAEC MPP 4','AE HITAEC MPP 6')

  all_prot_AE <- 0
    dat2_proteome <- Uniq_prot_proteome[,c(2,6:29)]
    dat2_proteome[dat2_proteome == 0] <- NA
    rownames(dat2_proteome) <- dat2_proteome[,1]
    dat2_proteome <- dat2_proteome[,-1]
    tdat_proteome <- t(dat2_proteome)
    dat_knn1_proteome <- impute.knn(tdat_proteome, k = 5)
    dat_knn_proteome <- t(dat_knn1_proteome$data)
    dat_log_proteome <- log2(dat_knn_proteome+1)
    all_prot_AE <- data.frame(normalizeQuantiles(dat_log_proteome))
    all_prot_AE$'Gene name' <- Uniq_prot_proteome[,1]
    all_prot_AE$'Protein' <- rownames(all_prot_AE)
    all_prot_AE <- all_prot_AE[,c(25,26,1:24)]
    rownames(all_prot_AE) <- rownames(Uniq_prot_proteome)
  
    colnames(all_prot_AE) <- c('Gene name',
                                      'Protein',
                                      'AE HCAEC CPP-P 2','AE HCAEC CPP-P 4','AE HCAEC CPP-P 6','AE HCAEC CPP-S 2','AE HCAEC CPP-S 4','AE HCAEC CPP-S 6','AE HCAEC DPBS 2','AE HCAEC DPBS 4','AE HCAEC DPBS 6','AE HCAEC MPP 2','AE HCAEC MPP 4','AE HCAEC MPP 6',
                                      'AE HITAEC CPP-P 2','AE HITAEC CPP-P 4','AE HITAEC CPP-P 6','AE HITAEC CPP-S 2','AE HITAEC CPP-S 4','AE HITAEC CPP-S 6','AE HITAEC DPBS 2','AE HITAEC DPBS 4','AE HITAEC DPBS 6','AE HITAEC MPP 2','AE HITAEC MPP 4','AE HITAEC MPP 6')
                                      
  all_AE_SandP <- merge(all_secr_AE, all_prot_AE, by = c('Gene name','Protein'))
  
  AE_heatmap <- all_AE_SandP[all_AE_SandP$Protein %in% prot_heatmap,]
  
  {
  means_AE_secr <- all_secr_AE[1:2]
    means_AE_secr$'MeanAE HCAEC CPP-P' <- rowMeans(all_secr_AE[,c(3:5)])
    means_AE_secr$'MeanAE HCAEC CPP-S' <- rowMeans(all_secr_AE[,c(6:8)])
    means_AE_secr[['MeanAE HCAEC (CPP-P + CPP-S)']] <- rowMeans(all_secr_AE[,c(3:8)])
    means_AE_secr$'MeanAE HCAEC DPBS' <- rowMeans(all_secr_AE[,c(9:10)])
    means_AE_secr$'MeanAE HCAEC MPP' <- rowMeans(all_secr_AE[,c(11:12)])
    means_AE_secr[['MeanAE HCAEC (DPBS + MPP)']] <- rowMeans(all_secr_AE[,c(9:12)])
    means_AE_secr$'MeanAE HITAEC CPP-P' <- rowMeans(all_secr_AE[,c(13:15)])
    means_AE_secr$'MeanAE HITAEC CPP-S' <- rowMeans(all_secr_AE[,c(16:18)])
    means_AE_secr[['MeanAE HITAEC (CPP-P + CPP-S)']] <-rowMeans(all_secr_AE[,c(13:18)])
    means_AE_secr$'MeanAE HITAEC DPBS' <- rowMeans(all_secr_AE[,c(19:21)])
    means_AE_secr$'MeanAE HITAEC MPP' <- rowMeans(all_secr_AE[,c(22:24)])
    means_AE_secr[['MeanAE HITAEC (DPBS + MPP)']] <-rowMeans(all_secr_AE[,c(19:24)])
  }
    
  {
  means_AE_prot <- all_prot_AE[1:2]
    means_AE_prot$'MeanAE HCAEC CPP-P' <- rowMeans(all_prot_AE[,c(3:5)])
    means_AE_prot$'MeanAE HCAEC CPP-S' <- rowMeans(all_prot_AE[,c(6:8)])
    means_AE_prot[['MeanAE HCAEC (CPP-P + CPP-S)']] <- rowMeans(all_prot_AE[,c(3:8)])
    means_AE_prot$'MeanAE HCAEC DPBS' <- rowMeans(all_prot_AE[,c(9:11)])
    means_AE_prot$'MeanAE HCAEC MPP' <- rowMeans(all_prot_AE[,c(12:14)])
    means_AE_prot[['MeanAE HCAEC (DPBS + MPP)']] <- rowMeans(all_prot_AE[,c(9:14)])
    means_AE_prot$'MeanAE HITAEC CPP-P' <- rowMeans(all_prot_AE[,c(15:17)])
    means_AE_prot$'MeanAE HITAEC CPP-S' <- rowMeans(all_prot_AE[,c(18:20)])
    means_AE_prot[['MeanAE HITAEC (CPP-P + CPP-S)']] <-rowMeans(all_prot_AE[,c(15:20)])
    means_AE_prot$'MeanAE HITAEC DPBS' <- rowMeans(all_prot_AE[,c(21:23)])
    means_AE_prot$'MeanAE HITAEC MPP' <- rowMeans(all_prot_AE[,c(24:26)])
    means_AE_prot[['MeanAE HITAEC (DPBS + MPP)']] <-rowMeans(all_prot_AE[,c(21:26)])
  }
   
  { 
  means_AE_SandP <- all_AE_SandP[,c(1:2)]
    means_AE_SandP$'Secr MeanAE HCAEC CPP-P' <- rowMeans(all_AE_SandP[,c(3:5)])
    means_AE_SandP$'Secr MeanAE HCAEC CPP-S' <- rowMeans(all_AE_SandP[,c(6:8)])
    means_AE_SandP$'Secr MeanAE HCAEC (CPP-P + CPP-S)' <- rowMeans(all_AE_SandP[,c(3:8)])
    means_AE_SandP$'Secr MeanAE HCAEC DPBS' <- rowMeans(all_AE_SandP[,c(9:10)])
    means_AE_SandP$'Secr MeanAE HCAEC MPP' <- rowMeans(all_AE_SandP[,c(11:12)])
    means_AE_SandP$'Secr MeanAE HCAEC (DPBS + MPP)' <- rowMeans(all_AE_SandP[,c(9:12)])
    means_AE_SandP$'Secr MeanAE HITAEC CPP-P' <- rowMeans(all_AE_SandP[,c(13:15)])
    means_AE_SandP$'Secr MeanAE HITAEC CPP-S' <- rowMeans(all_AE_SandP[,c(16:18)])
    means_AE_SandP$'Secr MeanAE HITAEC (CPP-P + CPP-S)' <- rowMeans(all_AE_SandP[,c(13:18)])
    means_AE_SandP$'Secr MeanAE HITAEC DPBS' <- rowMeans(all_AE_SandP[,c(19:21)])
    means_AE_SandP$'Secr MeanAE HITAEC MPP' <- rowMeans(all_AE_SandP[,c(22:24)])
    means_AE_SandP$'Secr MeanAE HITAEC (DPBS + MPP)' <- rowMeans(all_AE_SandP[,c(19:24)])
    means_AE_SandP$'Secr MeanAE' <- rowMeans(all_AE_SandP[,c(3:24)])
    means_AE_SandP$'Prot MeanAE HCAEC CPP-P' <- rowMeans(all_AE_SandP[,c(25:27)])
    means_AE_SandP$'Prot MeanAE HCAEC CPP-S' <- rowMeans(all_AE_SandP[,c(28:30)])
    means_AE_SandP$'Prot MeanAE HCAEC (CPP-P + CPP-S)' <- rowMeans(all_AE_SandP[,c(25:30)])
    means_AE_SandP$'Prot MeanAE HCAEC DPBS' <- rowMeans(all_AE_SandP[,c(31:33)])
    means_AE_SandP$'Prot MeanAE HCAEC MPP' <- rowMeans(all_AE_SandP[,c(34:36)])
    means_AE_SandP$'Prot MeanAE HCAEC (DPBS + MPP)' <- rowMeans(all_AE_SandP[,c(31:36)])
    means_AE_SandP$'Prot MeanAE HITAEC CPP-P' <- rowMeans(all_AE_SandP[,c(37:39)])
    means_AE_SandP$'Prot MeanAE HITAEC CPP-S' <- rowMeans(all_AE_SandP[,c(40:42)])
    means_AE_SandP$'Prot MeanAE HITAEC (CPP-P + CPP-S)' <- rowMeans(all_AE_SandP[,c(37:42)])
    means_AE_SandP$'Prot MeanAE HITAEC DPBS' <- rowMeans(all_AE_SandP[,c(43:45)])
    means_AE_SandP$'Prot MeanAE HITAEC MPP' <- rowMeans(all_AE_SandP[,c(46:48)])
    means_AE_SandP$'Prot MeanAE HITAEC (DPBS + MPP)' <- rowMeans(all_AE_SandP[,c(43:48)])
    means_AE_SandP$'Prot MeanAE' <- rowMeans(all_AE_SandP[,c(25:48)])
  }
    
  {
  means_AE_heatmap <- AE_heatmap[,c(1:2)]
    means_AE_heatmap$'Secr meanAE HCAEC (CPP-P + CPP-S)' <- rowMeans(AE_heatmap[,c(3:8)])
    means_AE_heatmap$'Secr meanAE HCAEC (DPBS + MPP)' <- rowMeans(AE_heatmap[,c(9:12)])
    means_AE_heatmap$'Secr meanAE HITAEC (CPP-P + CPP-S)' <- rowMeans(AE_heatmap[,c(13:18)])
    means_AE_heatmap$'Secr meanAE HITAEC (DPBS + MPP)' <- rowMeans(AE_heatmap[,c(19:24)])
    means_AE_heatmap$'Secr meanAE' <- rowMeans(AE_heatmap[,c(3:24)])
    means_AE_heatmap$'Prot meanAE HCAEC (CPP-P + CPP-S)' <- rowMeans(AE_heatmap[,c(25:30)])
    means_AE_heatmap$'Prot meanAE HCAEC (DPBS + MPP)' <- rowMeans(AE_heatmap[,c(31:36)])
    means_AE_heatmap$'Prot meanAE HITAEC (CPP-P + CPP-S)' <- rowMeans(AE_heatmap[,c(37:42)])
    means_AE_heatmap$'Prot meanAE HITAEC (DPBS + MPP)' <- rowMeans(AE_heatmap[,c(43:48)])
    means_AE_heatmap$'Prot meanAE' <- rowMeans(AE_heatmap[,c(43:48)])
  }

#heatmaps top50 of overexpression proteins in ctrl and experimental groups in HITAEC and HCAEC
  
    prot_top50_HCA_KM <- rownames(head(means_AE_secr[order(means_AE_secr$`MeanAE HCAEC (DPBS + MPP)`,decreasing = T),],n=50))
    length(prot_top50_HCA_KM)
    
    prot_top50_HCA_CI <- rownames(head(means_AE_secr[order(means_AE_secr$`MeanAE HCAEC (CPP-P + CPP-S)`,decreasing = T),],n=50))
    length(prot_top50_HCA_CI)
    
    prot_top50_HIT_KM <- rownames(head(means_AE_secr[order(means_AE_secr$`MeanAE HITAEC (DPBS + MPP)`,decreasing = T),],n=50))
    length(prot_top50_HIT_KM)
    
    prot_top50_HIT_CI <- rownames(head(means_AE_secr[order(means_AE_secr$`MeanAE HITAEC (CPP-P + CPP-S)`,decreasing = T),],n=50))
    length(prot_top50_HIT_CI)
    
pheatmap(dat_norm2[prot_top50_HCA_KM,c(1:4)],
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = F,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         breaks = seq(18.5, 24, length.out = 101),
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm2[prot_top50_HCA_CI,c(11:16)],
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = F,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         breaks = seq(18.5, 24, length.out = 101),
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm2[prot_top50_HIT_KM,c(5:10)],
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = F,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         breaks = seq(18.5, 24, length.out = 101),
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

pheatmap(dat_norm2[prot_top50_HIT_CI,c(17:22)],
         annotation_col = fact[,c(2,5)],
         cutree_cols = 2,
         cluster_cols = FALSE,
         cluster_rows = F,
         cellwidth = 12,
         cellheight = 10,
         border_color = "black",
         breaks = seq(18.5, 24, length.out = 101),
         color = colorRampPalette(c("#00bfff", '#005aeb', '#240935', '#800080','#cd00cd',"#ff19ff"))(100))

    
#write data in excel book (suppl datafile 17 - unique proteins)


  Datafile_17 <- createWorkbook()
  
  addWorksheet(Datafile_17, "HCAEC (DPBS + MPP)")
  writeData(Datafile_17, "HCAEC (DPBS + MPP)", Uniq_prot_secretome_ctrl_HCA)
  
  addWorksheet(Datafile_17, "HITAEC (DPBS + MPP)")
  writeData(Datafile_17, "HITAEC (DPBS + MPP)", Uniq_prot_secretome_ctrl_HIT)
  
  addWorksheet(Datafile_17, "(DPBS + MPP) - HCAEC vs HITAEC")
  writeData(Datafile_17, "(DPBS + MPP) - HCAEC vs HITAEC", Uniq_prot_secretome_HCA_ctrl_vs_HIT_ctrl)
  
  addWorksheet(Datafile_17, "HCAEC (CPP-P + CPP-S)")
  writeData(Datafile_17, "HCAEC (CPP-P + CPP-S)", Uniq_prot_secretome_exp_HCA)
  
  addWorksheet(Datafile_17, "HITAEC (CPP-P + CPP-S)")
  writeData(Datafile_17, "HITAEC (CPP-P + CPP-S)", Uniq_prot_secretome_exp_HIT)
  
  addWorksheet(Datafile_17, "(DPBS + MPP) - HITAEC vs HCAEC")
  writeData(Datafile_17, "(DPBS + MPP) - HITAEC vs HCAEC", Uniq_prot_secretome_HIT_ctrl_vs_HCA_ctrl)
  
  saveWorkbook(Datafile_17, file = "Supplementary Datafile 17 (Unique proteins).xlsx", overwrite = TRUE)


#write data in excel book (suppl datafile 15 - analysis of areas)

  Datafile_15 <- createWorkbook()
  
  addWorksheet(Datafile_15, "Raw") #all areas & counts for secretome
  writeData(Datafile_15, "Raw", Uniq_prot_secretome)
  
  addWorksheet(Datafile_15, "Analysis") #means for peak areas to secretome
  writeData(Datafile_15, "Analysis", Uniq_prot_areas_secretome)
  
  addWorksheet(Datafile_15, "Range") #46 of proteins, which were selected for heatmap
  writeData(Datafile_15, "Range", all_heatmap[,c(1:25)])
  
  addWorksheet(Datafile_15, "Range anal") #intergroup means on Range
  writeData(Datafile_15, "Range anal", all_heatmap_means[,c(1:8)])
  
  saveWorkbook(Datafile_15, file = "Supplementary Datafile 15 (secretome across the samples).xlsx", overwrite = TRUE)

#write data in excel FULL book
    
  Unique_proteins <- createWorkbook()
  
  #lists with areas
  addWorksheet(Unique_proteins, "All Secretome") #all areas & counts for secretome
  writeData(Unique_proteins, "All Secretome", Uniq_prot_secretome)
  
  addWorksheet(Unique_proteins, "All Proteome") #all areas & counts for proteome
  writeData(Unique_proteins, "All Proteome", Uniq_prot_proteome)
  
  addWorksheet(Unique_proteins, "All areas Secr∩Prot") #All areas for secretome was intersection with proteome
  writeData(Unique_proteins, "All areas Secr∩Prot", all_protein)
  
  addWorksheet(Unique_proteins, "Areas heatmap") #All areas for heatmap in secretome was intersection with proteome
  writeData(Unique_proteins, "Areas heatmap", all_heatmap)
  
  addWorksheet(Unique_proteins, "Means areas Secr") #all means for all areas in secretome
  writeData(Unique_proteins, "Means areas Secr", Uniq_prot_areas_secretome)
  
  addWorksheet(Unique_proteins, "Means areas Prot") #all means for all areas in proteome
  writeData(Unique_proteins, "Means areas Prot", Uniq_prot_areas_proteome)
  
  addWorksheet(Unique_proteins, "Means areas S∩P") #all means for all areas in secretome was intersection with proteome
  writeData(Unique_proteins, "Means areas S∩P", all_protein_means)
  
  addWorksheet(Unique_proteins, "Means areas heatmap") #All means areas for heatmap's proteins in secretome was intersection with proteome
  writeData(Unique_proteins, "Means areas heatmap", all_heatmap_means)
  
  #lists with AE
  addWorksheet(Unique_proteins, "All Secr AE") #all AEs for secretome
  writeData(Unique_proteins, "All Secr AE", all_secr_AE)
  
  addWorksheet(Unique_proteins, "All Prot AE") #all AEs for proteome
  writeData(Unique_proteins, "All Prot AE", all_prot_AE)
  
  addWorksheet(Unique_proteins, "All AE Secr∩Prot") #All AEs for secretome was intersection with proteome
  writeData(Unique_proteins, "All AE Secr∩Prot", all_AE_SandP)
  
  addWorksheet(Unique_proteins, "AE heatmap") #All AEs for heatmap in secretome was intersection with proteome
  writeData(Unique_proteins, "AE heatmap", AE_heatmap)
  
  addWorksheet(Unique_proteins, "Means AE Secr") #all means for all AEs in secretome
  writeData(Unique_proteins, "Means AE Secr", means_AE_secr)
  
  addWorksheet(Unique_proteins, "Means AE Prot") #all means for all AEs in proteome
  writeData(Unique_proteins, "Means AE Prot", means_AE_prot)
  
  addWorksheet(Unique_proteins, "Means AE S∩P") #all means for all AEs in secretome was intersection with proteome
  writeData(Unique_proteins, "Means AE S∩P", means_AE_SandP)
  
  addWorksheet(Unique_proteins, "Means AE heatmap") #all means for AEs for heatmap's proteins in secretome was intersection with proteome
  writeData(Unique_proteins, "Means AE heatmap", means_AE_heatmap)
  
  #lists with unique proteins
  addWorksheet(Unique_proteins, "Secr Exp Proteins HIT")
  writeData(Unique_proteins, "Secr Exp Proteins HIT", Uniq_prot_secretome_exp_HIT)
  
  addWorksheet(Unique_proteins, "Secr Ctrl Proteins HIT")
  writeData(Unique_proteins, "Secr Ctrl Proteins HIT", Uniq_prot_secretome_ctrl_HIT)
  
  addWorksheet(Unique_proteins, "Secr Exp Proteins HCA")
  writeData(Unique_proteins, "Secr Exp Proteins HCA", Uniq_prot_secretome_exp_HCA)
  
  addWorksheet(Unique_proteins, "Secr Ctrl Proteins HCA")
  writeData(Unique_proteins, "Secr Ctrl Proteins HCA", Uniq_prot_secretome_ctrl_HCA)
  
  addWorksheet(Unique_proteins, "Secr Ctrl HCA vs Ctrl HIT")
  writeData(Unique_proteins, "Secr Ctrl HCA vs Ctrl HIT", Uniq_prot_secretome_HCA_ctrl_vs_HIT_ctrl)
  
  addWorksheet(Unique_proteins, "Secr Ctrl HIT vs Ctrl HCA")
  writeData(Unique_proteins, "Secr Ctrl HIT vs Ctrl HCA", Uniq_prot_secretome_HIT_ctrl_vs_HCA_ctrl)
  
  saveWorkbook(Unique_proteins, file = "Unique_proteins_sec+prot.xlsx", overwrite = TRUE)
  
  