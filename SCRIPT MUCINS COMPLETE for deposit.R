#################################################
# PART 1 - CORRELATION MATRIX

#This script contains 3 parts: correlation, survival and subgroup analysis

setwd("D:/Jonckheere/2020/R")
# set working directory

#packages
if(!require("reshape2")) install.packages("reshape2") ;
library(reshape2)
if(!require("ggplot2")) install.packages("ggplot2") ;
library(ggplot2)
if(!require("Hmisc")) install.packages("Hmisc") ;
library(Hmisc)
if(!require("PerformanceAnalytics")) install.packages("PerformanceAnalytics") ;
library(PerformanceAnalytics)
if(!require("FactoMineR")) install.packages("FactoMineR") ;
library(FactoMineR)
if(!require("factoextra")) install.packages("factoextra") ;
library(factoextra)

#retrieve RSEM value of the gene of interest from cbioportal.org
GENE <-read.csv2("MUCIN PAAD COX RSEM.csv")
head(GENE)

#Select gene of interest. Remove genes without value
GENEb <- GENE[, c(2,4,6,7,8,9,10,11,12,13,15,17,18,19,21,22)]
head(GENEb)
summary(GENEb)

#Matrix
corGENE <- round(cor(GENEb), 2)
head(corGENE)

#library reshape2 
melted_corGENE <- melt(corGENE)
head(melted_corGENE)

library(ggplot2)
ggplot(data = melted_corGENE, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Lower triangle
get_lower_tri<-function(corGENE){
  corGENE[upper.tri(corGENE)] <- NA
  return(corGENE)
}
# Upper triangle
get_upper_tri <- function(corGENE){
  corGENE [lower.tri(corGENE)]<- NA
  return(corGENE)
}

upper_tri <- get_upper_tri(corGENE)
upper_tri

# Melt matrix
library(reshape2)
melted_corGENEup <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggheatmapGENE <- ggplot(data = melted_corGENEup, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
ggheatmapGENE

#Pearson on heatmap
ggheatmapGENE + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


#r and p values

library(Hmisc)
rcorr(as.matrix(GENEb), type=c("pearson","spearman"))

# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res<-rcorr(as.matrix(GENEb))
flattenCorrMatrix(round(res$r,2), round(res$P,5))

symnum(cor(GENEb), abbr.colnames=FALSE)


#hierarchical clustering

library("FactoMineR")
library("factoextra")



#continues data
library(FactoMineR)
# 1. ACP 
res.pca <- PCA(GENEb, ncp = 3, graph = T)
# 2. HCPC
res.hcpc <- HCPC(res.pca, graph = FALSE)

#dendrogramme
fviz_dend(res.hcpc, 
          cex = 0.4,                     # text size
          palette = "jco",               # color ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Rectangle
          rect_border = "jco",           # rectangle color
          labels_track_height = 0.8      # space for text
)
fviz_cluster(res.hcpc,
             repel = TRUE,            # avoid overlap
             show.clust.cent = TRUE, # cluster centers
             palette = "jco",         # colors, voir ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

head(res.hcpc$data.clust, 10)
#cluster characteristics
res.hcpc$desc.var$quanti
#cluster description:
res.hcpc$desc.axes$quanti
#representative individuals
res.hcpc$desc.ind$para

#get cluster list - adjust number of patients
head(res.hcpc$data.clust$clust, 168)

##################################################################
#PART 2 - SURVIVAL

library(survival)
library(survminer)
library(ggplot2)

#Create a csv file with cluster in a column and survival data (time and status)
#If necessary, filter cluster with less than 3 patients
SURVIE <-read.csv2("PAAD SURVIE CLUSTER2 MUCIN.csv")

fit <- survfit(Surv(time, status)~cluster,data=SURVIE)
view(fit)

#Log rank
survdiff(Surv(time,status)~cluster,data=SURVIE)

#test wilcoxon
survdiff(Surv(time,status)~cluster,data=SURVIE,rho=1)

plot(survfit(Surv(time, status)~cluster,data=SURVIE), main = "Survival curve")

KMsurv <- ggsurvplot(fit, data=SURVIE, main = "Survival curve", pval=TRUE,
                      font.main = 18,
                      font.x =  16,
                      font.y = 16,
                      font.tickslab = 14,
                      conf.int = TRUE,
                      risk.table = TRUE,
                      tables.height = 0.2,
                      tables.theme = theme_cleantable(),
                      palette = c("coral1", "#E7B800"))
KMsurv

coxmodel <- coxph(Surv(time,status) ~ cluster, data=SURVIE)
summary(coxmodel)


#################################################################"
# PART 3 - Analysis of clinical features 

if(!require("ggplot2")) install.packages("ggplot2") ;
library(ggplot2)
if(!require("dplyr")) install.packages("dplyr") ;
library(dplyr)

#Create a csv file with a cluster column and clinical features (Sex, diagnosis age, mutation count, stage, N, grade...)
PAAD <-read.csv2("PAAD CLINICAL 2 clusters.csv")
head(PAAD)
str(PAAD)


#transform as categorial value
PAAD$cluster <-as.factor(PAAD$cluster)
str(PAAD$cluster)

tab1 <- table(PAAD$cluster)
tab2 <- table(PAAD$Histo.Grade)
tab3 <- table(PAAD$Mutation.Count)
tab4 <- table(PAAD$Lymph.Node)
tab5 <- table(PAAD$Tumor.Stage)
tab6 <- table(PAAD$Sex)
tab7 <- table(PAAD$Diagnosis.Age)
prop.table(tab1)
prop.table(tab2)
prop.table(tab3)
prop.table(tab4)
prop.table(tab5)
prop.table(tab6)
prop.table(tab7)

#Tables of clinical features
tabHG <- table(PAAD$Histo.Grade,PAAD$cluster)
tabHG
prop.table(tabHG)

tabMC <- table(PAAD$Mutation.Count,PAAD$cluster)
tabMC
prop.table(tabMC)


tabLN <- table(PAAD$Lymph.Node,PAAD$cluster)
tabLN
prop.table(tabLN)


tabTS <- table(PAAD$Tumor.Stage,PAAD$cluster)
tabTS
prop.table(tabTS)

tabSex <- table(PAAD$Sex,PAAD$cluster)
tabSex
prop.table(tabSex)


chisq.test(tabHG)
chisq.test(tabTS)
chisq.test(tabSex)
chisq.test(tabLN)

#percentage plots
ggplot(PAAD, aes(x = cluster)) + geom_bar(aes(fill = Histo.Grade), position = "fill") +scale_y_continuous(labels = scales::percent_format())
ggplot(PAAD, aes(x = cluster)) + geom_bar(aes(fill = Tumor.Stage), position = "fill") +scale_y_continuous(labels = scales::percent_format())
ggplot(PAAD, aes(x = cluster)) + geom_bar(aes(fill = Lymph.Node), position = "fill")+scale_y_continuous(labels = scales::percent_format())
ggplot(PAAD, aes(x = cluster)) + geom_bar(aes(fill = Sex), position = "fill")+scale_y_continuous(labels = scales::percent_format())

#Boxplots
pMUT<-ggplot(PAAD, aes(x=cluster, y=Mutation.Count, fill=cluster)) +
  geom_boxplot(fill=c("red2", "goldenrod1"))
pMUT

pAGE<-ggplot(PAAD, aes(x=cluster, y=Diagnosis.Age, fill=cluster )) +
  geom_boxplot(fill=c("red2", "goldenrod1"))
pAGE

anova(lm(PAAD$Mutation.Count~PAAD$cluster))
anova(lm(PAAD$Diagnosis.Age~PAAD$cluster))
