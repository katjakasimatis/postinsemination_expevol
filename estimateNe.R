#load data
data.filtered <- read.table("kasimatis_genomicdata_filtered.txt", header = TRUE, sep = "\t")

data.filtered$line <- substr(data.filtered$sample, start = 1, stop = 6)

data.filtered$CHR[data.filtered$CHR=="I"]   <- "1"
data.filtered$CHR[data.filtered$CHR=="II"]  <- "2"
data.filtered$CHR[data.filtered$CHR=="III"] <- "3"
data.filtered$CHR[data.filtered$CHR=="IV"]  <- "4"
data.filtered$CHR[data.filtered$CHR=="V"]   <- "5"
data.filtered$CHR[data.filtered$CHR=="X"]   <- "6"


#create a dataframe for the population at generation 0 (i.e., ancestor) and keep only SNPs where both reference AND alternate alleles have coverage in nuclear genome
anc <- subset(data.filtered, line == "ANC_AA" & REF != 0 & ALT != 0 & CHR != "MtDNA", select = c("SNP", "CHR", "REF", "ALT", "coverage"))
names(anc) <- c("SNP", "CHR", "anc.ref", "anc.alt", "anc.coverage")


#create a dataframe for each of the replicate populations at generation 31 (i.e., evolved) and keep only SNPs where both reference AND alternate alleles have coverage in nuclear genome
A_3_CO <- subset(data.filtered, line == "A_3_CO" & REF != 0 & ALT !=0 & time==31 & CHR != "MtDNA", select = c("SNP", "REF", "ALT", "coverage"))
names(A_3_CO) <- c("SNP", "evo.ref", "evo.alt", "evo.coverage")

#merge ancestor and evolved replicate by SNP
merged <- merge(anc, A_3_CO, by.x = "SNP", by.y = "SNP")


#calculate Ne based on Waples 1989
#Plan II: the initial sample and individuals contributing to the next generation are independent binomial samples (i.e., remove individuals before they contribute to next generation)

estimateNeWII <- function(p0, pt, cov0, covt, gen, ploidy){
    
  A <- ((p0/cov0 - pt/covt)^2) / ((0.5 * (p0/cov0 + pt/covt)) - (p0/cov0 * pt/covt))
  
  B <- (1/length(cov0)) * sum(A - ((1/cov0) + (1/covt)))
  
  NeWII <- -gen / (ploidy * log(1 - B))
  return(NeWII)
}

estimateNeWII(p0 = merged$anc.ref, pt = merged$evo.ref, cov0 = merged$anc.coverage, covt = merged$evo.coverage, gen = 30, ploidy = 2)


by.chromosome = matrix(nrow = 24, ncol = 7);
for (i in 1:24){
  name <- c("A_1_LA", "A_2_LA", "A_3_LA", "B_1_LA", "B_2_LA", "B_3_LA", "A_1_SO", "A_2_SO", "A_3_SO", "B_1_SO", "B_2_SO", "B_3_SO", "A_1_CO", "A_2_CO", "A_3_CO", "B_1_CO", 
            "B_2_CO", "B_3_CO", "A_1_SC", "A_2_SC", "A_3_SC", "B_1_SC", "B_2_SC", "B_3_SC")
  
  A <- subset(data.filtered, line == name[i] & REF != 0 & ALT !=0 & time==31 & CHR != "MtDNA", select = c("SNP", "REF", "ALT", "coverage"))
  names(A) <- c("SNP", "evo.ref", "evo.alt", "evo.coverage")
  
  B <- merge(anc, A, by.x = "SNP", by.y = "SNP")
  
  for (j in 1:6){
    Ne <- estimateNeWII(p0 = subset(B, CHR==(j))$anc.ref, pt = subset(B, CHR==(j))$evo.ref, cov0 = subset(B, CHR==(j))$anc.coverage, covt = subset(B, CHR==(j))$evo.coverage, 
                        gen = 30, ploidy = 2)
    
    by.chromosome[i, j+1] <- Ne
  }
 
  by.chromosome[i, 1] <- name[i]
}

by.chromosome <- as.data.frame(by.chromosome)
names(by.chromosome) <- c("line", "I", "II", "III", "IV", "V", "X")
by.chromosome$I   <- as.numeric(by.chromosome$I)
by.chromosome$II  <- as.numeric(by.chromosome$II)
by.chromosome$III <- as.numeric(by.chromosome$III)
by.chromosome$IV  <- as.numeric(by.chromosome$IV)
by.chromosome$V   <- as.numeric(by.chromosome$V)
by.chromosome$X   <- as.numeric(by.chromosome$X)


#calculate genome-wide mean per line
by.chromosome$mean <- rowSums(by.chromosome[, 2:7])/6


#Analysis of variance: Is there a significant effect of treatment on Ne?
by.chromosome$regime <- substr(by.chromosome$line, start = 5, stop = 6)
model1 <- aov(mean ~ regime, data = by.chromosome)


#T-test: Is there a significant difference between autosomal Ne and sex chromosome Ne?
autosome <- c(by.chromosome[, 2], by.chromosome[, 3], by.chromosome[, 4], by.chromosome[, 5], by.chromosome[, 6])
xchrom   <- by.chromosome[, 7]

t.test(by.chromosome[, 7], by.chromosome[, 2:6], alternative = "two.sided")


#What is the upper bound on Ne assuming unequal sex ratio, such that all females mate but not all males?
Nm <- function(Ne, Nf) {(Nf * Ne) / (4 * Nf - Ne)}
#let Nf = 2500
#let Ne range from 16.3-23.5% of 5000 (815-1175)

