############################################################################
####  Code for Jesšovnik A, Blažević I, Lemić D, Pajač Živković I         ###
####             ANT FAUNA OF ANNUAL AND PERENNIAL CROPS                 ###
####                                2019                                 ###
############################################################################

setwd("~/Desktop/R-agr")

install.packages("psych")
install.packages("vegan")
install.packages("ape")

library(readr)
library(lattice)
library(picante)
library(psych)

# RICHNESS
# read in data with number of species per plot
dataset.div <- read_csv("div.csv")
dataset.div$group <- as.factor(dataset.div$group) 

# set an order of treatments from natural to most intensive
levels(dataset.div$group)
dataset.div$group <- ordered(dataset.div$group,
                             levels = c("control", "perennial", "annual"))
#descriptive statistics
summary(dataset.div)

#descriptive statistics by landuse
describeBy(dataset.div, dataset.div$group)

boxplot(spnum ~ group,
        data = dataset.div,
        boxwex=0.7,
        col=(c("darkslategrey","azure3","gold")),
        ylab="Species richness",
        xlab="Land use")

## ABUNDANCE
dataset.abund <- read_csv("abund.csv")
dataset.abund$group <- as.factor(dataset.abund$group)

levels(dataset.abund$group)
dataset.abund$group <- ordered(dataset.abund$group,
                               levels = c("control", "perennial", "annual"))
levels(dataset.abund$group)

summary(dataset.abund)
describeBy(dataset.abund, dataset.abund$group)

boxplot(abund ~ group,
        data = dataset.abund,
        boxwex=0.7,
        col=(c("darkslategrey","azure3","gold")),
        ylab="Ant abundance",
        xlab="Land use")


####  ANOVA ON ANT RICHNESS
dataset.logdiv <- dataset.div
dataset.logdiv$spnum <- log10(dataset.logdiv$spnum)
# dependant variable not normal therefore log transform the data 

dataset.logdiv$group <- as.factor(dataset.logdiv$group) 

levels(dataset.logdiv$group)
dataset.logdiv$group <- ordered(dataset.logdiv$group,
                         levels = c("control", "perennial", "annual"))

# Compute the analysis of variance
res.anova <- aov(spnum ~ group, data = dataset.logdiv)
# Summary of the analysis
summary(res.anova)

# post hoc 
TukeyHSD(res.anova)

# Check the anova assumptions: homogeneity of variance assumption
library(car)
leveneTest(spnum ~ group, data = dataset.logdiv)

#### ANOVA ON ANT ABUNDANCE

dataset.logabund <- dataset.abund
dataset.logabund$abund <- log10(dataset.logabund$abund)
View(dataset.logabund)

dataset.logabund$group <- as.factor(dataset.logabund$group) 

levels(dataset.logabund$group)
dataset.logabund$group <- ordered(dataset.logabund$group,
                                levels = c("control", "perennial", "annual"))

is.factor(dataset.logabund$group)

# Compute the analysis of variance
res.anova3 <- aov(abund ~ group, data = dataset.logabund)
# Summary of the analysis
summary(res.anova3)

leveneTest(spnum ~ group, data = dataset.logdiv)

# post hoc - which groups are different?
TukeyHSD(res.anova3)


###  COMMUNITY ANALYSES

# for data use plot IDs as rownames (first column of data) and 
# species names as colnames
dataset.comm <- read.csv("comm.csv", header = TRUE, row.names = 1)

# total abundance in each sample
apply(dataset.comm, 1, sum)

# attach the land use 
landuse <- read_csv("comm-meta.csv")
landuse$group <- as.factor(landuse$group)
attach(landuse)

levels(landuse$group)
landuse$group <- ordered(landuse$group,
                         levels = c("control", "perennial", "annual"))
      
# Shannon-Weiner Index
H <- diversity(dataset.comm, index = "shannon") 
summary(H) #gives summary statistics for the plots
View(H)

boxplot(H ~ landuse$group,
        data = dataset.comm,
        boxwex=0.7,
        col=(c("darkslategrey","azure3","gold")),
        ylab="Ant Diversity (H')",
        xlab="Land use")

# summary stats
divH <- read_csv("div-H.csv")
divH$group <- as.factor(divH$group)
divH$group <- ordered(divH$group,
                         levels = c("control", "perennial", "annual"))
describeBy(divH, divH$group)

# Shannon per land use

dataset.group <- read.csv("comm-groups.csv", header = TRUE, row.names = 1)
Hg <- diversity(dataset.group, index = "shannon") 
summary(Hg) 
View(Hg)

#### Other indices
D <- diversity(dataset.comm, index = "simpson") # Simpson Index
plot(D ~ landuse$group)

J <- H/log(specnumber(dataset.comm)) # Pielou evenness
plot(J ~ landuse$group)

jaccard <- vegdist(dataset.comm, method = "jaccard") # Jaccard distance
jaccard

## Fisher’s alpha
fish.a<-fisher.alpha(dataset.comm, MARGIN = 1)
fish.a
View(fish.a)

## saved as a csv file "div-f", for next step

#######################################
# ANOVA DIVERSITY

logdivH <- read_csv("div-f.csv")
logdivH$fisher <- log10(logdivH$fisher)

logdivH$group <- as.factor(logdivH$group) 
levels(logdivH$group)
logdivH$group <- ordered(logdivH$group,
                                 levels = c("control", "perennial", "annual"))
boxplot(fisher ~ group,
        data = logdivH,
        boxwex=0.7,
        col=(c("darkslategrey","azure3","gold")),
        ylab="log Fisher",
        xlab="Land use")

res.anova2 <- aov(fisher ~ group, data = logdivH)
summary(res.anova2)

leveneTest(spnum ~ group, data = logdivH)

TukeyHSD(res.anova2)

####### NMDS Non-metric Multidimensional Scaling

set.seed(3)
data.mds<-metaMDS(dataset.comm, distance = "bray", k = 2, trymax = 100, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE) 

stressplot(data.mds)

# make grey poligons that conect plots of same land use
treat=c(rep("control",3),rep("perennial",3), rep("annual",3))
ordiplot(data.mds,type="n")
ordihull(data.mds,groups=treat,draw="polygon",col="grey90",label=F)

# plot just the plots, colour by land use, pch=19 means plot a circle
points(data.mds, "sites", pch = 19, cex=1.5, col = "darkslategrey", select = landuse$group == 
         "control")
points(data.mds, "sites", pch = 19, cex=1.5, col = "azure3", select = landuse$group == 
         "perennial")
points(data.mds, "sites", pch = 19, cex=1.5, col = "gold", select = landuse$group == 
         "annual")

orditorp(data.mds,display="species",col="darkgrey",air=0.01)
orditorp(data.mds,display="sites",col=c(rep("darkslategrey",3),rep("azure3",3),rep("gold",3)),
         air=0.01,cex=1.5)





