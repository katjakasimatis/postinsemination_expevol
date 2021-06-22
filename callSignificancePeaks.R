library(dplyr)

#load time series GLM data
la <- read.table("FileS5_Model2_Sheet1.txt", header = TRUE, sep = "\t")
la.genomic <- subset(la, CHR != "MtDNA")
#la = within strain pre and post (WS-P&P)

co <- read.table("FileS5_Model2_Sheet2.txt", header = TRUE, sep = "\t")
co.genomic <- subset(co, CHR != "MtDNA")
#co = between strain pre and post (BS-P&P)

so <- read.table("FileS5_Model2_Sheet3.txt", header = TRUE, sep = "\t")
so.genomic <- subset(so, CHR != "MtDNA")
#so = within strain post only (WS-PO)

sc <- read.table("FileS5_Model2_Sheet4.txt", header = TRUE, sep = "\t")
sc.genomic <- subset(sc, CHR != "MtDNA")
#sc = between strain post only (BS-PO)


#Bonferronfi cut-off
la.bonf <- 1 / (0.5 * length(la.genomic$SNP))
co.bonf <- 1 / (0.5 * length(co.genomic$SNP))
so.bonf <- 1 / (0.5 * length(so.genomic$SNP))
sc.bonf <- 1 / (0.5 * length(sc.genomic$SNP)) 


#define significance peaks: at least 5 significant SNPs within a 1kb window
chrs <- c("I", "II", "III", "IV", "V", "X")

la.peaks <- c();
for (j in 1:length(chrs)) {
  
  A <- subset(la.genomic, CHR == chrs[j] & pvalue < la.bonf)
  
  for (i in 1:length(A$POS)) {
    A <- A[order(A$POS, decreasing = FALSE), ]
    
    if(A$SNP[i] %in% la.peaks$SNP == FALSE){
      B <- seq(A$POS[i], A$POS[i] + 999, by = 1)
      
      if(length(A[A$POS %in% B, 3]) >= 5){
        C <- A[A$POS %in% B, ]
        C$peak <- interaction(sprintf("chr_%s", j), sprintf("peak_%s", i), sep = "_")
        
        la.peaks <- rbind(la.peaks, C)
      }
      
      else{
        D <- A[A$POS %in% B, ]
        D$peak <- NA
        
        la.peaks <- rbind(la.peaks, D)
      }  
    }
    
    else{
      Z <- NA
    }
    
  }  
}
#number of significance peaks: 55
la.onlyp <- subset(la.peaks, is.na(peak)==FALSE)
la.summary <- c();
for (i in 1:length(levels(as.factor(la.onlyp$peak)))) {
  A <- levels(as.factor(la.onlyp$peak))[i]
  
  start <- min(la.onlyp[la.onlyp$peak %in% A, 2])
  end   <- max(la.onlyp[la.onlyp$peak %in% A, 2])
  name  <- interaction(A, "LA", sep = "_")
  
  la.summary <- rbind(la.summary, data.frame(CHR      = la.onlyp[la.onlyp$peak %in% A, 1][1],
                                             startPos = start,
                                             stopPos  = end,
                                             peak     = name))
}

co.peaks <- c();
for (j in 1:length(chrs)) {
  
  A <- subset(co.genomic, CHR == chrs[j] & pvalue < co.bonf)
  
  for (i in 1:length(A$POS)) {
    A <- A[order(A$POS, decreasing = FALSE), ]
    
    if(A$SNP[i] %in% co.peaks$SNP == FALSE){
      B <- seq(A$POS[i], A$POS[i] + 999, by = 1)
      
      if(length(A[A$POS %in% B, 3]) >= 5){
        C <- A[A$POS %in% B, ]
        C$peak <- interaction(sprintf("chr_%s", j), sprintf("peak_%s", i), sep = "_")
        
        co.peaks <- rbind(co.peaks, C)
      }
      
      else{
        D <- A[A$POS %in% B, ]
        D$peak <- NA
        
        co.peaks <- rbind(co.peaks, D)
      }  
    }
    
    else{
      Z <- NA
    }
    
  }  
}
#number of significance peaks: 81
co.onlyp <- subset(co.peaks, is.na(peak)==FALSE)
co.summary <- c();
for (i in 1:length(levels(as.factor(co.onlyp$peak)))) {
  A <- levels(as.factor(co.onlyp$peak))[i]
  
  start <- min(co.onlyp[co.onlyp$peak %in% A, 2])
  end   <- max(co.onlyp[co.onlyp$peak %in% A, 2])
  name  <- interaction(A, "CO", sep = "_")
  
  co.summary <- rbind(co.summary, data.frame(CHR      = co.onlyp[co.onlyp$peak %in% A, 1][1],
                                             startPos = start,
                                             stopPos  = end,
                                             peak     = name))
}

so.peaks <- c();
for (j in 1:length(chrs)) {
  
  A <- subset(so.genomic, CHR == chrs[j] & pvalue < so.bonf)
  
  for (i in 1:length(A$POS)) {
    A <- A[order(A$POS, decreasing = FALSE), ]
    
    if(A$SNP[i] %in% so.peaks$SNP == FALSE){
      B <- seq(A$POS[i], A$POS[i] + 999, by = 1)
      
      if(length(A[A$POS %in% B, 3]) >= 5){
        C <- A[A$POS %in% B, ]
        C$peak <- interaction(sprintf("chr_%s", j), sprintf("peak_%s", i), sep = "_")
        
        so.peaks <- rbind(so.peaks, C)
      }
      
      else{
        D <- A[A$POS %in% B, ]
        D$peak <- NA
        
        so.peaks <- rbind(so.peaks, D)
      }  
    }
    
    else{
      Z <- NA
    }
    
  }  
}
#number of significance peaks: 113
so.onlyp <- subset(so.peaks, is.na(peak)==FALSE)
so.summary <- c();
for (i in 1:length(levels(as.factor(so.onlyp$peak)))) {
  A <- levels(as.factor(so.onlyp$peak))[i]
  
  start <- min(so.onlyp[so.onlyp$peak %in% A, 2])
  end   <- max(so.onlyp[so.onlyp$peak %in% A, 2])
  name  <- interaction(A, "SO", sep = "_")
  
  so.summary <- rbind(so.summary, data.frame(CHR      = so.onlyp[so.onlyp$peak %in% A, 1][1],
                                             startPos = start,
                                             stopPos  = end,
                                             peak     = name))
}

sc.peaks <- c();
for (j in 1:length(chrs)) {
  
  A <- subset(sc.genomic, CHR == chrs[j] & pvalue < sc.bonf)
  
  for (i in 1:length(A$POS)) {
    A <- A[order(A$POS, decreasing = FALSE), ]
    
    if(A$SNP[i] %in% sc.peaks$SNP == FALSE){
      B <- seq(A$POS[i], A$POS[i] + 999, by = 1)
      
      if(length(A[A$POS %in% B, 3]) >= 5){
        C <- A[A$POS %in% B, ]
        C$peak <- interaction(sprintf("chr_%s", j), sprintf("peak_%s", i), sep = "_")
        
        sc.peaks <- rbind(sc.peaks, C)
      }
      
      else{
        D <- A[A$POS %in% B, ]
        D$peak <- NA
        
        sc.peaks <- rbind(sc.peaks, D)
      }  
    }
    
    else{
      Z <- NA
    }
    
  }  
}
#number of significance peaks: 118
sc.onlyp <- subset(sc.peaks, is.na(peak)==FALSE)
sc.summary <- c();
for (i in 1:length(levels(as.factor(sc.onlyp$peak)))) {
  A <- levels(as.factor(sc.onlyp$peak))[i]
  
  start <- min(sc.onlyp[sc.onlyp$peak %in% A, 2])
  end   <- max(sc.onlyp[sc.onlyp$peak %in% A, 2])
  name  <- interaction(A, "SC", sep = "_")
  
  sc.summary <- rbind(sc.summary, data.frame(CHR      = sc.onlyp[sc.onlyp$peak %in% A, 1][1],
                                             startPos = start,
                                             stopPos  = end,
                                             peak     = name))
}


#join peak files
peak.summary <- sc.summary %>% full_join(so.summary, by = "startPos") %>% full_join(co.summary, by = "startPos") %>% full_join(la.summary, by = "startPos")


