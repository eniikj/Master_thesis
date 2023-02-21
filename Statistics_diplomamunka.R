#1#
#setting the working directory
setwd("C:/Users/Eniko Juhasz/brca_tcga_files/rnaseq_files")



#2#
#reading in the tables needed
#brca norm seq data
norm_BRCA_data_Eniko <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/norm_BRCA_data_Eniko.txt", row.names = 1)
#coding mutation
Breast_kodolo_mutaciok <- read.delim("C:/Users/Eniko Juhasz/brca_tcga_files/my files/Breast_kodolo_mutaciok.txt", row.names=1)
#disruptive mutation
Breast_rombolo_mutaciok <- read.delim("C:/Users/Eniko Juhasz/brca_tcga_files/my files/Breast_rombolo_mutaciok.txt", row.names=1)


#3#
#extracting the data for the TP53 gene
Breast_kodolo_mutaciok<- data.frame(t(Breast_kodolo_mutaciok))
coding_TP53 <- Breast_kodolo_mutaciok[,"TP53", drop=FALSE]
coding_TP53$TCGA <- rownames(coding_TP53)

Breast_rombolo_mutaciok <- data.frame(t(Breast_rombolo_mutaciok))
disruptive_TP53 <- Breast_rombolo_mutaciok[,"TP53", drop=FALSE]
disruptive_TP53$TCGA <- rownames(disruptive_TP53)




#4#
#ratio of the coding and disruptive mutations regarding TP53
mutant_coding_TP53<- coding_TP53[coding_TP53$TP53 =="1", ]
wild_coding_TP53<- coding_TP53[coding_TP53$TP53 =="0", ]
mutant_coding_list <- list(mutant_coding_TP53$TCGA)
wild_coding_list <- list(wild_coding_TP53$TCGA)

mutant_disruptive_TP53<- disruptive_TP53[disruptive_TP53$TP53 =="1", ]
wild_disruptive_TP53<- disruptive_TP53[disruptive_TP53$TP53 =="0", ]
mutant_disruptive_list <- list(mutant_disruptive_TP53$TCGA)
wild_disruptive_list <- list(wild_disruptive_TP53$TCGA)


#5#
#coding expr tables
class(norm_BRCA_data_Eniko[1,1])

norm_BRCA_data_Eniko[,1:979]<-sapply(norm_BRCA_data_Eniko[,1:979],as.numeric)
library(tidyverse, lib.loc = "C:/Program Files/R/R-4.0.3/library")

brca_wild_coding <- select(norm_BRCA_data_Eniko,
                        -all_of(unlist(mutant_coding_list)))
write.csv(brca_wild_coding, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_wild_coding.txt")
brca_mutant_coding <- select(norm_BRCA_data_Eniko,
                             -all_of(unlist(wild_coding_list)))
write.csv(brca_mutant_coding, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_mutant_coding.txt")

#disruptive expr tables
brca_wild_disruptive <- select(norm_BRCA_data_Eniko,
                            -all_of(unlist(mutant_disruptive_list)))
write.csv(brca_wild_disruptive, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_wild_disruptive.txt")

brca_mutant_disruptive <- select(norm_BRCA_data_Eniko,
                             -all_of(unlist(wild_disruptive_list)))
write.csv(brca_mutant_disruptive, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_mutant_disruptive.txt")




#6#
#loading in the coding expression tables
brca_wild_coding <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_wild_coding.txt", row.names = 1)
brca_mutant_coding <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_mutant_coding.txt", row.names=1)


#binding the coding tables
brca_wild_coding <- rbind(c(0), brca_wild_coding)
rownames(brca_wild_coding)[rownames(brca_wild_coding) == "1"] <- "Mutation_status"
brca_wild_coding <- t(brca_wild_coding)


brca_mutant_coding <- rbind(c(1), brca_mutant_coding)
rownames(brca_mutant_coding)[rownames(brca_mutant_coding) == "1"] <- "Mutation_status"
brca_mutant_coding <- t(brca_mutant_coding)


coding <- rbind(brca_wild_coding, brca_mutant_coding)
coding <- t(coding)



#7#
#Mann Whitney U test
#CODING
coding2 <- coding[-1, ]
vector <- as.numeric(coding[1, ])
table(vector)

wilcox.test(as.numeric(coding2[1,]) ~ vector)
coding2[1, ]




#performing the test on every row for the coding mutations
out <- lapply(1:34315, function(x) wilcox.test(as.numeric(coding2[x, ]) ~ vector))
names(out) <- rownames(coding2)

#extracting the p values and putting them into a table
p.values <- lapply(out, `[`, c('p.value'))
p.values <- data.frame(matrix(unlist(p.values), nrow=length(p.values), byrow=T))

p.values <- data.frame(p.values)
rownames(p.values) <- rownames(coding2)
colnames(p.values) <- "p_value"


#calculating the averages for wild and for mutant
wild <- data.frame(rowMeans(subset(coding2, select = 1:674), na.rm = TRUE))
colnames(wild) <- "mean_wild"
mutant <- data.frame(rowMeans(subset(coding2, select = 675:979), na.rm = TRUE))
colnames(mutant) <- "mean_mutant"


#calculating the medians
wild_median <- data.frame(apply(coding2[ ,1:674], 1, median, na.rm = T))
colnames(wild_median) <- "median_wild"
mutant_median <- data.frame(apply(coding2[ ,675:979], 1, median, na.rm = T))
colnames(mutant_median) <- "median_mutant"


#merging the tables
coding_sum <- data.frame(p.values, wild, mutant, wild_median, mutant_median)

#calculating the mean fold
coding_sum$FC <- coding_sum$median_mutant/coding_sum$median_wild

#adding the average and filtering below 100 av. expression
average <- data.frame(rowMeans(coding2))
colnames(average)  <- "average"
coding_sum <- data.frame(coding_sum, average)


coding_sum <- subset(coding_sum, average > 100)

#kerek?t?sek
coding_sum2 <- names(coding_sum)[2:8]
coding_sum[,(coding_sum2)] <- round(coding_sum[,(coding_sum2)], digits = 2)



#renaming FC to original_FC and newFC to FC
names(coding_sum)[names(coding_sum) == "FC"] <- "original_FC"
names(coding_sum)[names(coding_sum) == "newFC"] <- "FC"

p_valuevectorc <- as.numeric(coding_sum$p_value)
adjustedp <- data.frame(p.adjust(p_valuevectorc, "fdr"))
coding_sum$p_value <- adjustedp$p.adjust.p_valuevectorc...fdr..

coding_sum3 <- names(coding_sum)[1]
coding_sum[,(coding_sum3)] <- data.frame(signif(coding_sum[,(coding_sum3)], digits=3))

write.csv(coding_sum, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/coding_sum.txt")
#loading in the coding table
coding_sum <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/coding_sum.txt", row.names = 1)




#8#
#loading in the disruptive expression tables
brca_wild_disruptive <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_wild_disruptive.txt", row.names = 1)
brca_mutant_disruptive <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/brca_mutant_disruptive.txt", row.names=1)


#binding the disruptive tables
brca_wild_disruptive <- rbind(c(0), brca_wild_disruptive)
rownames(brca_wild_disruptive)[rownames(brca_wild_disruptive) == "1"] <- "Mutation_status"
brca_wild_disruptive <- t(brca_wild_disruptive)


brca_mutant_disruptive <- rbind(c(1), brca_mutant_disruptive)
rownames(brca_mutant_disruptive)[rownames(brca_mutant_disruptive) == "1"] <- "Mutation_status"
brca_mutant_disruptive <- t(brca_mutant_disruptive)


disruptive <- rbind(brca_wild_disruptive, brca_mutant_disruptive)
disruptive <- t(disruptive)



#9#
#Mann Whitney U test
#DISRUPTIVE
disruptive2 <- disruptive[-1, ]
vector2 <- as.numeric(disruptive[1, ])
table(vector2)

wilcox.test(as.numeric(disruptive2[1, ]), vector2)
disruptive2[1, ]

#performing the test on every row for the coding mutations
out2 <- lapply(1:34315, function(x) wilcox.test(as.numeric(disruptive2[x, ]) ~ vector2))
names(out2) <- rownames(disruptive2)

#extracting the p values and putting them into a table
p.values2 <- lapply(out2, `[`, c('p.value'))
p.values2 <- data.frame(matrix(unlist(p.values2), nrow=length(p.values2), byrow=T))

p.values2 <- data.frame(p.values2)
rownames(p.values2) <- rownames(disruptive2)
colnames(p.values2) <- "p_value"

#calculating the averages
wild2 <- data.frame(rowMeans(subset(disruptive2, select = 1:875), na.rm = TRUE))
colnames(wild2) <- "mean_wild"
mutant2 <- data.frame(rowMeans(subset(disruptive2, select = 876:979)))
colnames(mutant2) <- "mean_mutant"


wild_median2 <- data.frame(apply(disruptive2[ ,1:875], 1, median, na.rm = T))
colnames(wild_median2) <- "median_wild"
mutant_median2 <- data.frame(apply(disruptive2[ ,876:979], 1, median, na.rm = T))
colnames(mutant_median2) <- "median_mutant"

#merging the tables
disruptive_sum <- data.frame(p.values2, wild2, mutant2, wild_median2, mutant_median2)

#calculating the mean fold
disruptive_sum$FC <- disruptive_sum$median_mutant/disruptive_sum$median_wild


#adding the average and filtering below 100 av. expression
average2 <- data.frame(rowMeans(disruptive2))
colnames(average2)  <- "average"
disruptive_sum <- data.frame(disruptive_sum, average2)


disruptive_sum <- subset(disruptive_sum, average > 100)


#kerek?t?sek
disruptive_sum2 <- names(disruptive_sum)[2:8]
disruptive_sum[,(disruptive_sum2)] <- round(disruptive_sum[,(disruptive_sum2)], digits = 2)




#renaming FC to original_FC and newFC to FC
names(disruptive_sum)[names(disruptive_sum) == "FC"] <- "original_FC"
names(disruptive_sum)[names(disruptive_sum) == "newFC"] <- "FC"


p_valuevectord <- as.numeric(disruptive_sum$p_value)
adjustedp2 <- data.frame(p.adjust(p_valuevectord, "fdr"))
disruptive_sum$p_value <- adjustedp2$p.adjust.p_valuevectord...fdr..

disruptive_sum3 <- names(disruptive_sum)[1]
disruptive_sum[,(disruptive_sum3)] <- signif(disruptive_sum[,(disruptive_sum3)], digits=3)


write.csv(disruptive_sum, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/disruptive_sum.txt")
#loading in the disruptive table
disruptive_sum <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/disruptive_sum.txt", row.names = 1)




#10#
#determining up-, and downregulating genes
#adding the reciprocal function
reciprocal <- function(x){
  y <- 1/x
  return(y)
}


#CODING#
coding_sum$newFC <- with(coding_sum, ifelse(mean_wild < mean_mutant, FC, reciprocal(FC)))
#adding new columns for up and downregulation
coding_sum$regulation <- with(coding_sum, ifelse(mean_wild < mean_mutant, "up", "down"))
#FC cutoff
coding_cutoff <- subset(coding_sum, newFC > 2)

#DISRUPTIVE#
disruptive_sum$newFC <- with(disruptive_sum, ifelse(mean_wild < mean_mutant, FC, reciprocal(FC)))
#adding new columns for up and downregulation
disruptive_sum$regulation <- with(disruptive_sum, ifelse(mean_wild < mean_mutant, "up", "down"))
#FC cutoff
disruptive_cutoff <- subset(disruptive_sum, newFC > 2)



#12#
#filtering based on p value <.01 and FC <... or >...
#CODING
coding_filter <- subset(coding_sum, p_value < 0.01)
coding_filter <- subset(coding_filter, FC > 2)


#DISRUPTIVE
disruptive_filter <- subset(disruptive_sum, p_value < 0.01)
disruptive_filter <- subset(disruptive_filter, FC > 2)

#filtering it for the coding genes
names <- rownames(coding_filter)
disruptive_filter <- subset(disruptive_filter, rownames(disruptive_filter) %in% names)



#saving and loading in coding_filter and disruptive_filter
write.csv(coding_filter, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/coding_filter.txt")
coding_filter <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/coding_filter.txt", row.names = 1)

coding_filter2 <- coding_filter[order(coding_filter$p_value),]

write.csv(disruptive_filter, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/disruptive_filter.txt")
disruptive_filter <- read.csv("C:/Users/Eniko Juhasz/brca_tcga_files/my files/disruptive_filter.txt", row.names = 1)




#13#
#gene ontology
library(dplyr, lib.loc = "C:/Program Files/R/R-4.0.3/library")
library(tidyverse)

#CODING UP
coding_up <- filter(coding_cutoff, 
                    regulation=="up")

cu_genes <- data.frame(rownames(coding_up))
write.xlsx(cu_genes, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/GO_codingdup.xlsx")

                          

coding_down <- filter(coding_cutoff,
                      regulation=="down")
cd_genes <- data.frame(rownames(coding_down))
write.xlsx(cd_genes, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/GO_codingdown.xlsx")

#DISRUPTIVE
disr_up <- filter(disruptive_cutoff,
                  regulation=="up" )
du_genes <- data.frame(rownames(disr_up))
write.xlsx(du_genes, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/GO_disrup.xlsx")


disr_down <- filter(disruptive_cutoff,
                  regulation=="down" )
dd_genes <- data.frame(rownames(disr_down))
write.xlsx(dd_genes, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/GO_disrdown.xlsx")



#14#
#CHOOSING 5 GENES FROM EACH
#CODING UPREGULATION, DOWNREGULATION
#ordering by pvalues
codingfilterup <- filter(coding_filter,
                         regulation == "up")
codingfilterdown <- filter(coding_filter,
                           regulation =="down")


codingfilterup <- codingfilterup[order(codingfilterup$p_value),]
genes_codingup <- codingfilterup[1:5,]


codingfilterdown <- codingfilterdown[order(codingfilterdown$p_value),]
genes_codingdown <- codingfilterdown[1:5,]

#DISRUPTIVE UPREGULATION, DOWNREGULATION
#ordering by pvalues
disruptivefilterup <- filter(disruptive_filter,
                         regulation == "up")
disruptivefilterdown <- filter(disruptive_filter,
                           regulation =="down")


disruptivefilterup <- disruptivefilterup[order(disruptivefilterup$p_value),]
genes_disrup <- disruptivefilterup[1:5, ]



disruptivefilterdown <- disruptivefilterdown[order(disruptivefilterdown$p_value),]
genes_disrdown <- disruptivefilterdown[1:5,]

codingfilterup2 <- data.frame(rownames(codingfilterup))
codingfilterdown2 <- data.frame(rownames(codingfilterdown))
disruptivefilterup2 <- data.frame(rownames(disruptivefilterup))
disruptivefilterdown2 <- data.frame(rownames(disruptivefilterdown))

write.xlsx(codingfilterup2, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/codingfilterup.xlsx")
write.xlsx(codingfilterdown2, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/codingfilterdown.xlsx")
write.xlsx(disruptivefilterup2, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/disrfilterup.xlsx")
write.xlsx(disruptivefilterdown2, "C:/Users/Eniko Juhasz/brca_tcga_files/my files/disrfilterdown.xlsx")







#15#
#Box plot
mutation_status <- data.frame(coding[1,])
names(mutation_status)[names(mutation_status) == "coding.1..."] <- "Mutation_status"
mutation_status$Mutation_status <- with(mutation_status, ifelse(Mutation_status=="0", "Vad t?pus", "Mut?ns t?pus"))

#ggplot2 package
library("ggplot2")

#CODING UPREGULATION
#1
#PSAT1

PSAT1 <- data.frame(t(subset(coding, rownames(coding) == "PSAT1")))
PSAT1 <- merge(PSAT1, mutation_status, by = 0, all = TRUE)
names(PSAT1)[names(PSAT1) == "PSAT1"] <- "dataPSAT1"


pPSAT1 <- ggplot(PSAT1, aes(x=as.factor(Mutation_status),
                            y=dataPSAT1,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "PSAT1 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.198e-55")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataPSAT1,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,10000)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pPSAT1  

#2
#CBX2
CBX2 <- data.frame(t(subset(coding, rownames(coding) == "CBX2")))
CBX2 <- merge(CBX2, mutation_status, by = 0, all = TRUE)

names(CBX2)[names(CBX2) == "CBX2"] <- "dataCBX2"

pCBX2 <- ggplot(CBX2, aes(x=as.factor(Mutation_status),
                            y=dataCBX2,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "CBX2 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 2.775e-54") +
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataCBX2,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,10000)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pCBX2

#3
#CDC20
CDC20 <- data.frame(t(subset(coding, rownames(coding) == "CDC20")))
CDC20 <- merge(CDC20, mutation_status, by = 0, all = TRUE)

names(CDC20)[names(CDC20) == "CDC20"] <- "dataCDC20"

pCDC20 <- ggplot(CDC20, aes(x=as.factor(Mutation_status),
                            y=dataCDC20,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "CDC20 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 1.305e-53") +
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataCDC20,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,20000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5))+
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pCDC20 

#4
#MYBL2
MYBL2 <- data.frame(t(subset(coding, rownames(coding) == "MYBL2")))
MYBL2<- merge(MYBL2, mutation_status, by = 0, all = TRUE)

names(MYBL2)[names(MYBL2) == "MYBL2"] <- "dataMYBL2"

pMYBL2 <- ggplot(MYBL2, aes(x=as.factor(Mutation_status),
                            y=dataMYBL2,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "MYBL2 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 4.872e-52") +
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataMYBL2,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,22500)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pMYBL2 


#5
#CENPA
CENPA <- data.frame(t(subset(coding, rownames(coding) == "CENPA")))
CENPA <- merge(CENPA, mutation_status, by = 0, all = TRUE)

names(CENPA)[names(CENPA) == "CENPA"] <- "dataCENPA"

pCENPA <- ggplot(CENPA, aes(x=as.factor(Mutation_status),
                            y=dataCENPA,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "CENPA g?nexpresszi?ja",
       subtitle = "p ?rt?k = 5.166e-51") +
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataCENPA,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,3000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pCENPA 




#CODING DOWNREGULATION

#6
#MLPH
MLPH<- data.frame(t(subset(coding, rownames(coding) == "MLPH")))
MLPH <- merge(MLPH, mutation_status, by = 0, all = TRUE)
names(MLPH)[names(MLPH) == "MLPH"] <- "dataMLPH"


pMLPH <- ggplot(MLPH, aes(x=as.factor(Mutation_status),
                            y=dataMLPH,
                            fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "MLPH g?nexpresszi?ja",
       subtitle = "p ?rt?k = 4.138e-58")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataMLPH,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,90000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pMLPH  


#7
#TBC1D9
TBC1D9<- data.frame(t(subset(coding, rownames(coding) == "TBC1D9")))
TBC1D9 <- merge(TBC1D9, mutation_status, by = 0, all = TRUE)
names(TBC1D9)[names(TBC1D9) == "TBC1D9"] <- "dataTBC1D9"


pTBC1D9 <- ggplot(TBC1D9, aes(x=as.factor(Mutation_status),
                          y=dataTBC1D9,
                          fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "TBC1D9 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 1.582e-56")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataTBC1D9,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7) +
  ylim(0,125000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pTBC1D9  


#8
#THSD4

THSD4<- data.frame(t(subset(coding, rownames(coding) == "THSD4")))
THSD4 <- merge(THSD4, mutation_status, by = 0, all = TRUE)
names(THSD4)[names(THSD4) == "THSD4"] <- "dataTHSD4"


pTHSD4 <- ggplot(THSD4, aes(x=as.factor(Mutation_status),
                          y=dataTHSD4,
                          fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "THSD4 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.620e-55")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataTHSD4,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,40000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pTHSD4 


#9
#LINC00504

LINC00504<- data.frame(t(subset(coding, rownames(coding) == "LINC00504")))
LINC00504 <- merge(LINC00504, mutation_status, by = 0, all = TRUE)
names(LINC00504)[names(LINC00504) == "LINC00504"] <- "dataLINC00504"


pLINC00504 <- ggplot(LINC00504, aes(x=as.factor(Mutation_status),
                          y=dataLINC00504,
                          fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "LINC00504 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.403e-54")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataLINC00504,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,5000)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pLINC00504  


#10
#SRARP

SRARP<- data.frame(t(subset(coding, rownames(coding) == "SRARP")))
SRARP <- merge(SRARP, mutation_status, by = 0, all = TRUE)
names(SRARP)[names(SRARP) == "SRARP"] <- "dataSRARP"


pSRARP <- ggplot(SRARP, aes(x=as.factor(Mutation_status),
                          y=dataSRARP,
                          fill = Mutation_status)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "SRARP g?nexpresszi?ja",
       subtitle = "p ?rt?k = 1.641e-53")+
  geom_jitter(aes(x = as.factor(Mutation_status),
                  y = dataSRARP,
                  color = Mutation_status),
              size = 1,
              alpha = 0.7)+
  ylim(0,10000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pSRARP  

#DISRUPTIVE

mutation_status2 <- data.frame(disruptive[1,])
names(mutation_status2)[names(mutation_status2) == "disruptive.1..."] <- "Mutation_status2"
mutation_status2$Mutation_status2 <- with(mutation_status2, ifelse(Mutation_status2=="0", "Vad t?pus", "Mut?ns t?pus"))


#DISR UPREGULATION

#11
#MCM10

MCM10<- data.frame(t(subset(disruptive, rownames(disruptive) == "MCM10")))
MCM10 <- merge(MCM10, mutation_status2, by = 0, all = TRUE)
names(MCM10)[names(MCM10) == "MCM10"] <- "dataMCM10"


pMCM10 <- ggplot(MCM10, aes(x=as.factor(Mutation_status2),
                            y=dataMCM10,
                            fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "MCM10 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 5.011e-22")+
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataMCM10,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,5000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pMCM10


#12
#MYBL2

MYBL2<- data.frame(t(subset(disruptive, rownames(disruptive) == "MYBL2")))
MYBL2 <- merge(MYBL2, mutation_status2, by = 0, all = TRUE)
names(MYBL2)[names(MYBL2) == "MYBL2"] <- "dataMYBL2"


pMYBL2 <- ggplot(MYBL2, aes(x=as.factor(Mutation_status2),
                            y=dataMYBL2,
                            fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "MYBL2 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 1.499e-21")+
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataMYBL2,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,22000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pMYBL2


#13
#CENPA

CENPA2 <- data.frame(t(subset(disruptive, rownames(disruptive) == "CENPA")))
CENPA2 <- merge(CENPA2, mutation_status2, by = 0, all = TRUE)

names(CENPA2)[names(CENPA2) == "CENPA"] <- "dataCENPA2"

pCENPA2 <- ggplot(CENPA2, aes(x=as.factor(Mutation_status2),
                            y=dataCENPA2,
                            fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "CENPA g?nexpresszi?ja",
       subtitle = "p ?rt?k = 2.378e-21") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataCENPA2,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,2500) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pCENPA2

#14
#ZNF695

ZNF695 <- data.frame(t(subset(disruptive, rownames(disruptive) == "ZNF695")))
ZNF695 <- merge(ZNF695, mutation_status2, by = 0, all = TRUE)

names(ZNF695)[names(ZNF695) == "ZNF695"] <- "dataZNF695"

pZNF695 <- ggplot(ZNF695, aes(x=as.factor(Mutation_status2),
                              y=dataZNF695,
                              fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "ZNF695 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.173e-21") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataZNF695,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,1000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pZNF695

#15
#AURKB

AURKB <- data.frame(t(subset(disruptive, rownames(disruptive) == "AURKB")))
AURKB <- merge(AURKB, mutation_status2, by = 0, all = TRUE)

names(AURKB)[names(AURKB) == "AURKB"] <- "dataAURKB"

pAURKB <- ggplot(AURKB, aes(x=as.factor(Mutation_status2),
                              y=dataAURKB,
                              fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "AURKB g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.956e-21") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataAURKB,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,4500) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pAURKB

#DISR DOWNREGULATION

#16
#CAPN8

CAPN8 <- data.frame(t(subset(disruptive, rownames(disruptive) == "CAPN8")))
CAPN8 <- merge(CAPN8, mutation_status2, by = 0, all = TRUE)

names(CAPN8)[names(CAPN8) == "CAPN8"] <- "dataCAPN8"

pCAPN8 <- ggplot(CAPN8, aes(x=as.factor(Mutation_status2),
                            y=dataCAPN8,
                            fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "CAPN8 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 2.366e-37") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataCAPN8,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,18000)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pCAPN8

#17
#MLPH

MLPH2 <- data.frame(t(subset(disruptive, rownames(disruptive) == "MLPH")))
MLPH2 <- merge(MLPH2, mutation_status2, by = 0, all = TRUE)

names(MLPH2)[names(MLPH2) == "MLPH"] <- "dataMLPH2"

pMLPH2 <- ggplot(MLPH2, aes(x=as.factor(Mutation_status2),
                          y=dataMLPH2,
                          fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "MLPH g?nexpresszi?ja",
       subtitle = "p ?rt?k = 1.893e-24") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataMLPH2,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,90000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pMLPH2

#18
#NOSTRIN

NOSTRIN <- data.frame(t(subset(disruptive, rownames(disruptive) == "NOSTRIN")))
NOSTRIN <- merge(NOSTRIN, mutation_status2, by = 0, all = TRUE)

names(NOSTRIN)[names(NOSTRIN) == "NOSTRIN"] <- "dataNOSTRIN"

pNOSTRIN <- ggplot(NOSTRIN, aes(x=as.factor(Mutation_status2),
                            y=dataNOSTRIN,
                            fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "NOSTRIN g?nexpresszi?ja",
       subtitle = "p ?rt?k = 2.317e-24") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataNOSTRIN,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,8000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pNOSTRIN


#19
#SRARP

SRARP2 <- data.frame(t(subset(disruptive, rownames(disruptive) == "SRARP")))
SRARP2 <- merge(SRARP2, mutation_status2, by = 0, all = TRUE)

names(SRARP2)[names(SRARP2) == "SRARP"] <- "dataSRARP2"

pSRARP2 <- ggplot(SRARP2, aes(x=as.factor(Mutation_status2),
                                y=dataSRARP2,
                                fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "SRARP g?nexpresszi?ja",
       subtitle = "p ?rt?k = 3.138e-23") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataSRARP2,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,10000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pSRARP2


#20
#GATA3

GATA3 <- data.frame(t(subset(disruptive, rownames(disruptive) == "GATA3")))
GATA3 <- merge(GATA3, mutation_status2, by = 0, all = TRUE)

names(GATA3)[names(GATA3) == "GATA3"] <- "dataGATA3"

pGATA3 <- ggplot(GATA3, aes(x=as.factor(Mutation_status2),
                              y=dataGATA3,
                              fill = Mutation_status2)) + 
  geom_boxplot(outlier.color = NULL,
               outlier.shape = NA,
               alpha = 0.3) + 
  labs(x = "TP53 mut?ci?s st?tusza",
       y = "GATA3 g?nexpresszi?ja",
       subtitle = "p ?rt?k = 2.113e-22") +
  geom_jitter(aes(x = as.factor(Mutation_status2),
                  y = dataGATA3,
                  color = Mutation_status2),
              size = 1,
              alpha = 0.7)+
  ylim(0,200000) +
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_color_brewer(palette = "Dark2", direction = -1)+
  theme_minimal()+
  theme(legend.position="none",
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = c("Vad t?pus", "Mut?ns t?pus"))

pGATA3










#15#
#volcano plot
#CODING
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager',lib="C:/Program Files/R/R-4.0.3/library")

BiocManager::install('EnhancedVolcano',lib="C:/Program Files/R/R-4.0.3/library")
devtools::install_github('kevinblighe/EnhancedVolcano',lib="C:/Program Files/R/R-4.0.3/library")

library(EnhancedVolcano)

volcano_coding <- EnhancedVolcano(coding_plot,
                                  lab = rownames(coding_plot),
                                  x = 'log2Foldchange',
                                  y = 'pvalue',
                                  xlim = c(-6,6),
                                  pCutoff = 0.01,
                                  FCcutoff = log2(2),
                                  title = "Differenci?lisan expressz?lt g?nek a TP53 k?dol? mut?ci?ja eset?n",
                                  legendLabels=c('Nem szignigik?ns esetek','Fold Change alapj?n szignifik?ns esetek','P ?rt?k alapj?n szignifik?ns esetek',
                                                 'Szignifik?ns esetek'),
                                  legendPosition = 'right')

volcano_coding

coding_plot <- data.frame(log2Foldchange, pvalue)
rownames(coding_plot) = rownames(coding_sum)

log2Foldchange <- as.numeric(coding_sum$original_FC)
log2Foldchange = log2(log2Foldchange)
log2Foldchange <- as.numeric(log2Foldchange)
head(log2Foldchange)



pvalue <- as.numeric(coding_sum$p_value)
head(pvalue)

log2(1.44)



#DISRUPTIVE
log2Foldchange2 <- as.numeric(disruptive_sum$original_FC)
log2Foldchange2 = log2(log2Foldchange2)


pvalue2 <- as.numeric(disruptive_sum$p_value)

disruptive_plot <- data.frame(log2Foldchange2, pvalue2)
rownames(disruptive_plot) = rownames(disruptive_sum)

max(disruptive_plot$log2Foldchange2, na.rm = T)
min(-log10(disruptive_plot$pvalue), na.rm = T)

volcano_disr <- EnhancedVolcano(disruptive_plot,
                                lab = rownames(disruptive_plot),
                                x = 'log2Foldchange2',
                                y = 'pvalue2',
                                xlim = c(-6,6),
                                pCutoff = 0.01,
                                FCcutoff = log2(2),
                                title = "Differenci?lisan expressz?lt g?nek a TP53 diszrupt?v mut?ci?ja eset?n",
                                legendLabels=c('Nem szignigik?ns esetek','Fold Change alapj?n szignifik?ns esetek','P ?rt?k alapj?n szignifik?ns esetek',
                                               'Szignifik?ns esetek'),
                                legendPosition = 'right')

volcano_disr 
