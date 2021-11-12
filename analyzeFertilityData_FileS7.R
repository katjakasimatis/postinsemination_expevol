library(dplyr)
library(Hmisc)
library(lmerTest)
library(multcomp)
library(RColorBrewer)
library(ggplot2)
library(cowplot)

data <- read.csv("FileS7_FertilityData.csv", header = TRUE)

data$regime <- substr(data$line, start = 5, stop = 6)
data$regime[data$regime==""] <- "anc"

newdata <- as.data.frame(summarise(group_by(data, line, replicate, treatment), sum(total), sum(rfp.minus)))
names(newdata) <- c("line", "replicate", "treatment", "sum_total", "sum_rfp.minus")
newdata <- subset(newdata, replicate!="rep4")

newdata$selected_percent <- newdata$sum_rfp.minus / newdata$sum_total

### ANCESTRAL COMPETITIVE ABILITY ###
#Compute a 95% confidence interval and do a hypothesis test
#total reproductive success
ngm.counts <- subset(data, line=="Anc" & treatment=="NGM" & replicate != "rep4")
number.ngm.minus <- sum(ngm.counts$rfp.minus)
number.ngm.total <- sum(ngm.counts$total)

prop.test(x=number.ngm.minus, n=number.ngm.total, p=0.5, alternative="two.sided", conf.level = 0.95)  

#post-insemination reproductive success
aux.counts <- subset(data, line=="Anc" & treatment=="aux" & replicate != "rep4")
number.aux.minus <- sum(aux.counts$rfp.minus)
number.aux.total <- sum(aux.counts$total)

prop.test(x=number.aux.minus, n=number.aux.total, p=0.5, alternative="two.sided", conf.level = 0.95) 

anc.aux <- subset(newdata, treatment == "aux" & line == "Anc")
anc.aux.mean <- mean(anc.aux$selected_percent)

anc.ngm <- subset(newdata, treatment == "NGM" & line == "Anc")
anc.ngm.mean <- mean(anc.ngm$selected_percent)
####


### EVOLVED COMPETITIVE ABILITY ###
#total reproductive success: anc-evo
anc.evo.total <- glmer(cbind(rfp.plus, rfp.minus) ~ regime + (1|line), family = binomial(link = "logit"), data = subset(data, treatment == "NGM"))
anc.evo.mcp1 <- glht(anc.evo.total, linfct = mcp(regime = c("anc - LA = 0", "anc - CO = 0", "anc - SO = 0", "anc - SC = 0")))
summary(anc.evo.mcp1)


#post-insemination reproductive success: anc-evolved
anc.evo.sperm <- glmer(cbind(rfp.plus, rfp.minus) ~ regime + (1|line), family = binomial(link = "logit"), data = subset(data, treatment == "aux"))
anc.evo.mcp2 <- glht(anc.evo.sperm, linfct = mcp(regime = c("anc - LA = 0", "anc - CO = 0", "anc - SO = 0", "anc - SC = 0")))
summary(anc.evo.mcp2)


#total reproductive success: Is there an effects of experimental evolution under direct sexual selection?
ngm <- subset(data, treatment=="NGM" & line!="Anc")
total.glm <- glmer(cbind(rfp.plus, rfp.minus) ~ regime + (1|line), family = binomial(link = "logit"), data = ngm)
summary(total.glm)

#Type III Analysis of Variance Table
total.aov <- anova(total.glm)

#planned comparisons test
total.mcp <- glht(total.glm, linfct = mcp(regime = c("LA - SO = 0", "LA - CO = 0", "LA - SC = 0")))
summary(total.mcp)


#post-insemination reproductive success: Is there an effects of experimental evolution under direct sexual selection?
aux <- subset(data, treatment=="aux" & line!="Anc")
sperm.glm <- glmer(cbind(rfp.plus, rfp.minus) ~ regime + (1|line), family = binomial(link = "logit"), data = aux)
summary(sperm.glm)

#Type III Analysis of Variance Table
sperm.aov <- anova(sperm.glm)

#planned comparisons test
sperm.mcp <- glht(sperm.glm, linfct = mcp(regime = c("LA - SO = 0", "LA - CO = 0", "LA - SC = 0")))
summary(sperm.mcp)

