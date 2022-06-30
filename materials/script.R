# For this course, we follwed this material: file:///home/jcalle/Desktop/Statistical_Methods_for_Microbiome_Analysis/day1/rmarkdown_microbiome_with_results.html#1_Elements_of_a_phyloseq_object

default.par <- par()      
#par(default.par) # recover default graphical settings  

rm(list=ls())

#-------------------------------------#
#             LOAD LIBRARIES          #
#-------------------------------------#
library(phyloseq)
library(vegan)


#---------------------------#
#           DAY 1           #
#---------------------------#


#-------------------------------------#
#    LOAD, EXPLORE & PROCESS DATA     #
#-------------------------------------#
data(GlobalPatterns)
data <- GlobalPatterns

## nOTU x nSample
dim(data@otu_table)

## proportion of zeros
sum(data@otu_table==0)/(ntaxa(data)*nsamples(data))

## explore phyloseq object
data
head(data@otu_table)
View(data@sam_data)
head(data@tax_table)
data@phy_tree

## define SampleTypes of human origin for posterior analysis
human.samples <- c("Feces", "Skin", "Tongue")
data@sam_data$human <- rep(0,26)
data@sam_data[data@sam_data$SampleType %in% human.samples,"human"] <- 1
sample_data(data)$human<-factor(sample_data(data)$human,levels=c(0,1),labels=c("No","Yes"))

## define metadata
metadata <- as.data.frame(data@sam_data[,c(1,6,7,8)])

## explore abundance per samples (hist, boxplot of logs...)
sample_sums(data)
hist(sample_sums(data), breaks = 10)
boxplot(log(sample_sums(data))~metadata$SampleType)

## filter OTU by abundance
data<-prune_taxa(taxa_sums(data)>10, data)

## filter sample by abundance
data<-prune_samples(sample_sums(data)>10, data)
dim(data@otu_table)

## transform abundances (TSS)
propdata<-transform_sample_counts(data, function(x) x/sum(x))
propdata@otu_table[1:5,1:5]
sample_sums(propdata@otu_table) # check success of transformation

## agglomerate taxa at kingdom level (sum counts for each sample at taxa level?)
dataKingdom<-tax_glom(data, "Kingdom")
dataKingdom
taxa_names(dataKingdom)<-as.vector(dataKingdom@tax_table[,1])
taxa_names(dataKingdom)
rank_names(data)

propdataKingdom<-tax_glom(propdata, rank_names(data)[1])
propdataKingdom

dataPhylum<-tax_glom(data, rank_names(data)[2])
dataPhylum

propdataPhylum<-tax_glom(propdata, rank_names(data)[2])
propdataPhylum
taxa_names(propdataPhylum)<-tax_table(propdataPhylum)[,2]
?tax_glom

summary(taxa_sums(propdataPhylum))

## merge less abundant taxa information into one "artificial" taxa
less_freq_taxa<-which(taxa_sums(propdataPhylum)<0.8)
propdataPhylum2<-merge_taxa(propdataPhylum,less_freq_taxa)

## rename less frequent taxa as "rare"
taxa_names(propdataPhylum2)[taxa_names(propdataPhylum2)%in%taxa_names(propdataPhylum)[less_freq_taxa]]<-"rare"


#-------------------------------------#
#           ABUNDANCE PLOTS           #
#-------------------------------------#
## barplot at Kindom lvl
# plot_bar(dataKingdom, fill="Kingdom")
# plot_bar(propdataKingdom, fill="Kingdom")

## barplot at Phylum lvl
# plot_bar(propdataPhylum, fill="Phylum")

## barplot of agglomerated table with rare taxa
otu2<-otu_table(propdataPhylum2)[c(1:4,6,5),]

# par(xpd=TRUE, mar= c(5, 4, 4, 7) + 0.1, pin=c(6,4))
# barplot(otu2, ylim=c(0,1.1), col=2:7, beside=FALSE, space=0, las=3,
#         xlab="", ylab="%", main="Composition at Phylum level",xaxt="n")
# legend(26,1,legend=rownames(otu2), fill=2:7, bty="n")

par(default.par)

## boxplot for human samples - Phylum
# boxplot(as.data.frame(t(otu2)[sample_data(data)$human=="Yes",]), col=rainbow(6))
# title("Philum composition, human samples")
# 
# ## boxplot for non-human samples - PHylum
# boxplot(as.data.frame(t(otu2)[sample_data(data)$human=="No",]), col=rainbow(6))
# title("Philum composition, non human samples")

## heatmap
# plot_heatmap(dataPhylum)
# plot_heatmap(dataPhylum,sample.label="SampleType",  low="#66CCFF",high="#000033", na.value="white")
# 
# plot_heatmap(propdataPhylum2)
# plot_heatmap(propdataPhylum2,sample.label="SampleType",  low="#66CCFF",high="#000033", na.value="white")


#-------------------------------------#
#         PHYLOGENETIC TREE           #
#-------------------------------------#

## aglomerate at Class lvl
dataClass<-tax_glom(data, "Class")

## consider only a certain subset of taxa
d<-subset_taxa(dataClass, Phylum =="Acidobacteria")

## Phylogenetic tree for d
# plot_tree(d,ladderize="left",size="abundance",color="human",label.tips="Class") 



#---------------------------#
#           DAY 2           #
#---------------------------#

#-----------------------------------------------#
#                 ALPHA-DIVERSITY               #
#-----------------------------------------------#

## plot richness
rich <- estimate_richness(data, measures = c("Observed", "Chao1", "Shannon"))
rich
exp(rich$Shannon)[-c(1:15)] # efective value of species in the sample, easier to interpret than Shannon index

plot_richness(data, x="SampleType", measures=c("Observed","Chao1","Shannon"), color="human")

## shannon index between human and non-human samples
boxplot(rich$Shannon~metadata$human)

## the saturation of the number of species per sample size indicates that our coverage is adequate for our samples
rarecurve(t(otu_table(data))[1:5,],step=5000)

## richness acording to the data type
boxplot(rich$Shannon~metadata$SampleType)
title("Richness: Shannon by SampleType")

## pairwise comparison of shannon diversity index between sample types
pairwise.wilcox.test(rich$Shannon, metadata$SampleType, p.adjust.method ="fdr")
kruskal.test(rich$Shannon, metadata$SampleType)


#-----------------------------------------------#
#                 BETA-DIVERSITY                #
#-----------------------------------------------#

## distance matrices
  ## bray-curtis
  BC.dist<-distance(data, method ="bray")
  as.matrix(BC.dist)[1,2]
  ## unifrac
  UF.dist<-distance(data, method = "unifrac")
  
  
## subset only human data and calculate bray-curtis distance between samples
data.human<-subset_samples(data, human=="Yes")
BC.dist.human<-distance(data.human, method ="bray")

## perform MDS ordination of human samples with bray-curtis disntance and plot it
BC.MDS.human<-ordinate(data.human, "MDS", distance = BC.dist.human)
plot_ordination(data.human, BC.MDS.human,   color="SampleType")


#-----------------------------------------------#
#  MULTIVARIATE DIFFERENTIAL ABUNDANCE TESTING  #
#-----------------------------------------------#

## perform a PERMANOVA test via "adonis" function
test.adonis<-adonis2(BC.dist.human~sample_data(data.human)$SampleType)
test.adonis


#-----------------------------------------------#
#               coda4microbiome                 #
#-----------------------------------------------#

library("coda4microbiome")

#------------------------------------#
#           Binary Variable          #
#------------------------------------#

## seed to reproduce tutorial results
set.seed(123) 

## load the data
data(Crohn, package = "coda4microbiome")   

## check distribution of our binary variable (dependents variable)
table(y_Crohn)

## train model
coda_glmnet_Crohn<-coda_glmnet(x=x_Crohn,y=y_Crohn)

## number of taxa kept as relevant  
coda_glmnet_Crohn$taxa.num

## taxa coefficients
coda_glmnet_Crohn$`log-contrast coefficients`

## predictions (as the linear combination of obtained coefficients an log of relevant taxa)
coda_glmnet_Crohn$predictions

## perform a permutation test to assess overfitting on our study
null_acc<-coda_glmnet_null(x=x_Crohn,y=y_Crohn, niter=10)

## check accuracy of prediction in null model
null_acc$accuracy # for AUC ~0.5 it means there is no discrimintation, which is expected in a null model

## check confidence interval
null_acc$`confidence interval` # our model results do not overlap with the null modell, which means that our results are significant


#------------------------------------#
#        Continuous Variable         #
#------------------------------------#

## load continuous data
data(sCD14, package = "coda4microbiome")

## train the model
coda_glmnet_SCD14<-coda_glmnet(x=x_sCD14,y=y_sCD14, lambda = "lambda.min") # we take "lambda.min" because with lambda.1se we don't have ANY log-ratio

## the square of the correlation between (self) predicted and observed values
coda_glmnet_SCD14$`apparent Rsq`


#------------------------------------#
#             Log-Ratios             #
#------------------------------------#

Crohn_logratios<-explore_logratios(x=x_Crohn,y=y_Crohn, measure = "glm")






#--------------------------------------#

rm(list=ls())


#------------------------------------#
#       Longitudinal signature       #
#------------------------------------#

set.seed(123) # to reproduce the results

data(ecam_filtered, package = "coda4microbiome")   # load the data

x=x_ecam # microbiome abundance
x_time = metadata$day_of_life;    # observation times
subject_id = metadata$studyid;   # subject id
y= metadata$diet;           # diet ("bd"= breast diet, "fd"=formula diet)
ini_time = 0;
end_time = 90;
ecam_logratios<-explore_lr_longitudinal(x,y, x_time, subject_id, ini_time, end_time)

ecam_logratios$`order of importance`
ecam_logratios$`name of most important variables`
ecam_logratios$`max log-ratio`
taxanames[as.numeric(ecam_logratios$`max log-ratio`)]

plot_signature_curves(c(30,35), c(-1,1),x,y, x_time, subject_id, ini_time, end_time)


#------------------------------------#
#         Variable-selection         #
#------------------------------------#

set.seed(123)

ecam <-coda_glmnet_longitudinal (x,y, x_time, subject_id, ini_time, end_time, lambda="lambda.min",nfolds=4)
