#upload libraries
library(plyr)
library(valr)
library(dplyr)
library(tidyverse)
library(multcomp)
library(emmeans)


#load the raw data
data <- read.csv("Documents/UO.academics/Phillips_Lab/experiments/experimental_evolution/2018_FullExperiment/genomic_data/2020.05.25_expevol_genomics_merged_allelicdepth.txt", header = TRUE, sep = "\t")

#remove repetitive regions
bed <- read.table("Desktop/repeats.bed", header = FALSE)
names(bed) <- c("chrom", "start", "end")
bed <- as_tibble(bed)

data$chrom <- data$CHROM
data$start <- data$POS - 1
data$end   <- data$POS
data2 <- as_tibble(data)

data.masked <- as.data.frame(bed_subtract(data2, bed))
data.masked <- data.masked[, 1:77]

expevol <- c();
for(i in 5:ncol(data.masked)){
  expevol <- rbind(expevol, data.frame(CHR = as.character(data.masked$CHROM),
                                       POS = as.numeric(as.character(data.masked$POS)),
                                       REFALLELE = as.character(data.masked$REF),
                                       ALTALLELE = as.character(data.masked$ALT),
                                       ALL = as.character(data.masked[,i]),
                                       sample = colnames(data.masked)[i]));
}

expevol <- expevol[!grepl("/", expevol$ALL), ]
expevol <- expevol[expevol$REFALLELE %in% c("A","T","G","C") & expevol$ALTALLELE %in% c("A","T","G","C"), ]

expevol$REF <- gsub(",[0-9]+", "", perl=T, expevol$ALL)
expevol$ALT <- gsub("[0-9]+,", "", perl=T, expevol$ALL)
expevol$SNP <- paste(expevol$CHR, expevol$POS, expevol$REFALLELE, expevol$ALTALLELE, sep="_")

expevol$regime <- substr(expevol$sample, start = 5, stop = 6)
#LA = within strain pre- and post (WS-P&P)
#SO = within strain post only (WS-PO)
#CO = between strain pre- and post (BS-P&P)
#SC = between strain post only (BS-PO)

expevol$time <- substr(expevol$sample, start = 8, stop = 9)
expevol$time <- as.character(expevol$time)
expevol$time[expevol$time=="s0"]   <- "0"
expevol$time[expevol$time=="s4"]   <- "13"
expevol$time[expevol$time=="s7"]   <- "22"
expevol$time[expevol$time=="s1"]   <- "31"
expevol$time <- as.numeric(expevol$time)
#0 = generation 0 (i.e., ancestor)
#13 = generation 13 (i.e., after selective event 4)
#s7 = generation 22 (i.e., after selective event 7)
#s10 = generation 31 (i.e., evolved)


#estimate the coverage
expevol$REF <- as.numeric(as.character(expevol$REF))
expevol$ALT <- as.numeric(as.character(expevol$ALT))
expevol$coverage <- expevol$REF + expevol$ALT


#estimate observed nucleotide diversity in the ancestral population (following Begun et al. 2007. PLoS Biology)
ancestor <- subset(expevol, regime == "AA" & coverage >=20 & CHR != "MtDNA")

anc.pi <- c();
for (i in 1:length(ancestor$SNP)){
  print(i)
  A <- ancestor[i, ];
  
  cov <- A$coverage
  alt <- A$ALT
  
  B <- (2 * alt *(cov - alt)) / (cov * (cov - 1))
  
  C <- data.frame(CHR = A$CHR,
                  POS = A$POS,
                  SNP = A$SNP,
                  pi  = B)
  anc.pi <- rbind(anc.pi, C)
}

polymorphic <- subset(anc.pi, pi != 0)


#remove repetitive regions
bed <- read.table("repeats.bed", header = FALSE)
names(bed) <- c("chrom", "start", "end")
bed <- subset(bed, chrom != "MtDNA")
bed <- as_tibble(bed)

anc.pi$chrom <- anc.pi$CHR
anc.pi$start <- anc.pi$POS - 1
anc.pi$end   <- anc.pi$POS
anc.pi2 <- as_tibble(anc.pi)

anc.pi.masked <- as.data.frame(bed_subtract(anc.pi2, bed))
anc.pi.masked <- anc.pi.masked[, 1:4]

polymorphic.masked <- subset(anc.pi.masked, pi != 0)


#calculate Watterson's theta
#define 1Kb windows for each chromosome
chrI  <- seq(0, 15069000, by = 1000)  #(length 15,069 windows)
chrII <- seq(0, 15274000, by = 1000)  #(length 15,274 windows)
chrIII <- seq(0, 13779000, by = 1000)  #(length 13,790 windows)
chrIV <- seq(0, 17491000, by = 1000)  #(length 17,491 windows)
chrV  <- seq(0, 20916000, by = 1000)  #(length 20,916 windows)
chrX  <- seq(0, 17717000, by = 1000)  #(length 17,717 windows)


#get number of SNPs per window
slide_theta <- c();
for (i in 2:length(chrI)){
  #print(i);
  
  A <- subset(ancestor, POS >= chrI[i] - 10000 & POS < chrI[i] & CHR == "I")
  out <- data.frame(CHR    = A$CHR[1],
                    domain = A$domain[1],
                    nSNPs  = length(A$SNP))
  
  slide_theta <- rbind(slide_theta, out)
}
slide_theta <- subset(slide_theta, is.na(CHR) == FALSE)


#calculate the harmonic number (denominator value)
harmonicNumber <- 0
poolsize       <- 2500 * 3
numChromosomes <- poolsize * 2
for (i in 1:(numChromosomes - 1)) {
  harmonicNumber = harmonicNumber + 1.0/i
}
#for pool size of 7,500 worms: harmonic number = 10.19299

slide_theta$theta <- slide_theta$nSNP / harmonicNumber
slide_theta$theta_perBP <- slide_theta$theta / 1000
slide_theta$bp <- c(1:nrow(slide_theta))


ks.test(subset(slide_theta, domain != "center" & CHR == "X")$theta_perBP, subset(slide_theta, domain == "center" & CHR == "X")$theta_perBP)


#estimate the SFS
ancestor$minor <- pmin(ancestor$REF, ancestor$ALT)
ancestor$MAF <- ancestor$minor / ancestor$coverage

ancestor$AAF <- ancestor$ALT / ancestor$coverage

ancestor$domain <- c(rep("left.arm", length(subset(ancestor, CHR == "I" & POS <= 3858000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "I" & POS > 3858000 & POS <= 11040000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "I" & POS > 11040000)$SNP)),
                     rep("left.arm", length(subset(ancestor, CHR == "II" & POS <= 4879000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "II" & POS > 4879000 & POS <= 12020000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "II" & POS > 12020000)$SNP)),
                     rep("left.arm", length(subset(ancestor, CHR == "III" & POS <= 3722000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "III" & POS > 3722000 & POS <= 10340000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "III" & POS > 10340000)$SNP)),
                     rep("left.arm", length(subset(ancestor, CHR == "IV" & POS <= 3896000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "IV" & POS > 3896000 & POS <= 12970000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "IV" & POS > 12970000)$SNP)),
                     rep("left.arm", length(subset(ancestor, CHR == "V" & POS <= 5897000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "V" & POS > 5897000 & POS <= 16550000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "V" & POS > 16550000)$SNP)),
                     rep("left.arm", length(subset(ancestor, CHR == "X" & POS <= 6137000)$SNP)),
                     rep("center", length(subset(ancestor, CHR == "X" & POS > 6137000 & POS <= 12480000)$SNP)),
                     rep("right.arm", length(subset(ancestor, CHR == "X" & POS > 12480000)$SNP)))


ks.test(subset(ancestor, domain != "center" & CHR == 1)$pi, subset(ancestor, domain == "center" & CHR == 1)$pi)
