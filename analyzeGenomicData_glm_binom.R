#upload libraries
library(plyr)
library(valr)
library(dplyr)
library(tidyverse)
library(multcomp)
library(emmeans)


#load the raw data
data <- read.csv("kasimatis_genomics_merged_allelicdepth.txt", header = TRUE, sep = "\t")

#remove repetitive regions
bed <- read.table("repeats.bed", header = FALSE)
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

#filter the data: keep SNPs with at least 20x coverage but less than 281x coverage (i.e., remove top 5%)
data.filtered  <- expevol[expevol$coverage >= 20, ]
data.filtered  <- data.filtered[data.filtered$coverage < 281, ]
data.filtered  <- data.filtered[complete.cases(data.filtered[ ,c(7,8)]), ]


####
#Model 1
#subset data for ancestor-evolved comparison
anc.evo  <- subset(data.filtered, time == 0 | time == 31)
factors2 <- count(as.factor(anc.evo$SNP))
anc.evo  <- anc.evo[anc.evo$SNP %in% as.character(factors2[factors2$freq >= 10,]$x), ]

anc.evo$regime <- as.factor(anc.evo$regime)

#use the GML with one term and binomial distribution
#planned comparison between ancestor and evolved
start.time <- Sys.time()
GLM_AncEvo_table <- c();
MCP_AncEvo_table <- c();
for (i in 1:length(levels(as.factor(anc.evo$SNP)))){
  print(i);
  
  #keep only SNPs with data for ancestor AND evolved lines
  A <- anc.evo[anc.evo$SNP == levels(as.factor(anc.evo$SNP))[i], ];
  
  prod <- prod(A$time)
  
  if (prod == 0){
    #run glm
    B <- glm(cbind(REF, ALT) ~ regime, family = binomial(link = "logit"), data = A)
    C <- summary(B)
    
    #check for errors in GML caused by missing data
    D <- try(c(C$coefficients[2,1], C$coefficients[2,4], C$coefficients[3,1], C$coefficients[3,4], C$coefficients[4,1], C$coefficients[4,4], C$coefficients[5,1], C$coefficients[5,4]))
    
    if (class(D) == "try-error") {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        regimeCO_slope = NA,
                        regimeCO_pval  = NA,
                        regimeLA_slope = NA,
                        regimeLA_pval  = NA,
                        regimeSC_slope = NA,
                        regimeSC_pval  = NA,
                        regimeSO_slope = NA,
                        regimeSO_pval  = NA)
    }
    else {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        regimeCO_slope = C$coefficients[2,1],
                        regimeCO_pval  = C$coefficients[2,4],
                        regimeLA_slope = C$coefficients[3,1],
                        regimeLA_pval  = C$coefficients[3,4],
                        regimeSC_slope = C$coefficients[4,1],
                        regimeSC_pval  = C$coefficients[4,4],
                        regimeSO_slope = C$coefficients[5,1],
                        regimeSO_pval  = C$coefficients[5,4])
    }
    GLM_AncEvo_table <- rbind(GLM_AncEvo_table, out)
    
    #planned comparisons test: ancestor versus evolved
    E <- try(summary(glht(B, linfct = mcp(regime = c("AA - CO = 0", "AA - LA = 0", "AA - SC = 0", "AA - SO = 0")))))
    
    if (class(E) == "try-error") {
      cmp <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        AA_CO_slope = NA,
                        AA_CO_pval  = NA,
                        AA_LA_slope = NA,
                        AA_LA_pval  = NA,
                        AA_SC_slope = NA,
                        AA_SC_pval  = NA,
                        AA_SO_slope = NA,
                        AA_SO_pval  = NA,
                        row.names = NULL)
    }
    else {
      cmp <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        AA_CO_slope = E$test$coefficients[1],
                        AA_CO_pval  = E$test$pvalues[1],
                        AA_LA_slope = E$test$coefficients[2],
                        AA_LA_pval  = E$test$pvalues[2],
                        AA_SC_slope = E$test$coefficients[3],
                        AA_SC_pval  = E$test$pvalues[3],
                        AA_SO_slope = E$test$coefficients[4],
                        AA_SO_pval  = E$test$pvalues[4],
                        row.names = NULL)
    }
    MCP_AncEvo_table <- rbind(MCP_AncEvo_table, cmp)
  }
  end.time <- Sys.time()
}

####


###
#Model 2
#subset data into regimes
la.data.filtered <- subset(data.filtered, regime == "AA" | regime == "LA")
la.factors       <- count(as.factor(la.data.filtered$SNP))
la.data.filtered <- la.data.filtered[la.data.filtered$SNP %in% as.character(la.factors[la.factors$freq >= 9,]$x), ]
#LA regime masked repetitive regions after filtering: 202,933 (62.1% of before filtering total)

co.data.filtered <- subset(data.filtered, regime == "AA" | regime == "CO")
co.factors       <- count(as.factor(co.data.filtered$SNP))
co.data.filtered <- co.data.filtered[co.data.filtered$SNP %in% as.character(co.factors[co.factors$freq >= 9,]$x), ]
#CO regime masked repetitive regions after filtering: 222,744 (68.2% of before filtering total)

so.data.filtered <- subset(data.filtered, regime == "AA" | regime == "SO")
so.factors       <- count(as.factor(so.data.filtered$SNP))
so.data.filtered <- so.data.filtered[so.data.filtered$SNP %in% as.character(so.factors[so.factors$freq >= 9,]$x), ]
#SO regime masked repetitive regions after filtering: 200,329 (61.3% of before filtering total)

sc.data.filtered <- subset(data.filtered, regime == "AA" | regime == "SC")
sc.factors       <- count(as.factor(sc.data.filtered$SNP))
sc.data.filtered <- sc.data.filtered[sc.data.filtered$SNP %in% as.character(sc.factors[sc.factors$freq >= 9,]$x), ]
#SC regime masked repetitive regions after filtering: 204,954 (62.7% of before filtering total)

#use GML for each regime separately with time as the variable and a binomial distribution
#note: uses full ancestor for each regime
GLM_LA_table <- c();
GLM_SO_table <- c();
GLM_CO_table <- c();
GLM_SC_table <- c();

for (i in 1:length(levels(as.factor(la.data.filtered$SNP)))) {
  print(i);
  
  #create per SNP subset of the data
  A <- la.data.filtered[la.data.filtered$SNP == levels(as.factor(la.data.filtered$SNP))[i], ];
  
  t <- try(A[A$regime == "AA", ]$regime <- "AA", silent = TRUE)
  
  if (class(t) == "try-error"){
    fail <- A$SNP[1]
  }
  else {
    
    #run glm LA
    model1 <- glm(cbind(REF, ALT) ~ time, family = binomial(link = "logit"), data = A)
    B <- summary(model1)
    
    #check for errors in GML caused by missing data
    C <- try(c(B$coefficients[2,1], B$coefficients[2,2], B$coefficients[2,3], B$coefficients[2,4], silent = TRUE))
    
    if (class(C) == "try-error") {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = NA,
                        estimate  = NA,
                        std_error = NA,
                        z_value   = NA,
                        pvalue    = NA)
    }
    else {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = B$coefficients[1,1],
                        estimate  = B$coefficients[2,1],
                        std_error = B$coefficients[2,2],
                        z_value   = B$coefficients[2,3],
                        pvalue    = B$coefficients[2,4])
    }
    GLM_LA_table <- rbind(GLM_LA_table, out)
  }    
}

for (i in 1:length(levels(as.factor(co.data.filtered$SNP)))) {
  print(i);
  
  #create per SNP subset of the data
  A <- co.data.filtered[co.data.filtered$SNP == levels(as.factor(co.data.filtered$SNP))[i], ];
  
  t <- try(A[A$regime == "AA", ]$regime <- "AA", silent = TRUE)
  
  if (class(t) == "try-error"){
    fail <- A$SNP[1]
  }
  else {
    
    #run glm CO
    model1 <- glm(cbind(REF, ALT) ~ time, family = binomial(link = "logit"), data = A)
    B <- summary(model1)
    
    #check for errors in GML caused by missing data
    C <- try(c(B$coefficients[2,1], B$coefficients[2,2], B$coefficients[2,3], B$coefficients[2,4], silent = TRUE))
    
    if (class(C) == "try-error") {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = NA,
                        estimate  = NA,
                        std_error = NA,
                        z_value   = NA,
                        pvalue    = NA)
    }
    else {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = B$coefficients[1,1],
                        estimate  = B$coefficients[2,1],
                        std_error = B$coefficients[2,2],
                        z_value   = B$coefficients[2,3],
                        pvalue    = B$coefficients[2,4])
    }
    GLM_CO_table <- rbind(GLM_CO_table, out)
  }    
}

for (i in 1:length(levels(as.factor(so.data.filtered$SNP)))) {
  print(i);
  
  #create per SNP subset of the data
  A <- so.data.filtered[so.data.filtered$SNP == levels(as.factor(so.data.filtered$SNP))[i], ];
  
  t <- try(A[A$regime == "AA", ]$regime <- "AA", silent = TRUE)
  
  if (class(t) == "try-error"){
    fail <- A$SNP[1]
  }
  else {
    
    #run glm SO
    model1 <- glm(cbind(REF, ALT) ~ time, family = binomial(link = "logit"), data = A)
    B <- summary(model1)
    
    #check for errors in GML caused by missing data
    C <- try(c(B$coefficients[2,1], B$coefficients[2,2], B$coefficients[2,3], B$coefficients[2,4], silent = TRUE))
    
    if (class(C) == "try-error") {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = NA,
                        estimate  = NA,
                        std_error = NA,
                        z_value   = NA,
                        pvalue    = NA)
    }
    else {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = B$coefficients[1,1],
                        estimate  = B$coefficients[2,1],
                        std_error = B$coefficients[2,2],
                        z_value   = B$coefficients[2,3],
                        pvalue    = B$coefficients[2,4])
    }
    GLM_SO_table <- rbind(GLM_SO_table, out)
  }    
}

for (i in 1:length(levels(as.factor(sc.data.filtered$SNP)))) {
  print(i);
  
  #create per SNP subset of the data
  A <- sc.data.filtered[sc.data.filtered$SNP == levels(as.factor(sc.data.filtered$SNP))[i], ];
  
  t <- try(A[A$regime == "AA", ]$regime <- "AA", silent = TRUE)
  
  if (class(t) == "try-error"){
    fail <- A$SNP[1]
  }
  else {
    
    #run glm SC
    model1 <- glm(cbind(REF, ALT) ~ time, family = binomial(link = "logit"), data = A)
    B <- summary(model1)
    
    #check for errors in GML caused by missing data
    C <- try(c(B$coefficients[2,1], B$coefficients[2,2], B$coefficients[2,3], B$coefficients[2,4], silent = TRUE))
    
    if (class(C) == "try-error") {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = NA,
                        estimate  = NA,
                        std_error = NA,
                        z_value   = NA,
                        pvalue    = NA)
    }
    else {
      out <- data.frame(CHR = A$CHR[1],
                        POS = A$POS[1],
                        SNP = A$SNP[1],
                        intercept = B$coefficients[1,1],
                        estimate  = B$coefficients[2,1],
                        std_error = B$coefficients[2,2],
                        z_value   = B$coefficients[2,3],
                        pvalue    = B$coefficients[2,4])
    }
    GLM_SC_table <- rbind(GLM_SC_table, out)
  }    
}

###
