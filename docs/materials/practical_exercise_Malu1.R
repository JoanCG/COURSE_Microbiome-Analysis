library(phyloseq)
library(vegan)
library(coda4microbiome)


# Data set Qin et al. (2012) Type2 Diabetes microbiome study:  https://www.nature.com/articles/nature11450
# 1.	Upload the data as a phyloseq object as follows: 
 data<-readRDS("qin_all.rds")
# 
# 2.	Use phyloseq library to preprocess and describe the dataset: 
#   a.	Number of taxa, number of samples, number of variables in metadata
 data
ntaxa(data) 
nsamples(data)
dim(data@sam_data)

 
# b.	Agglomerate at the genus level and normalize the data to relative abundances

rank_names(data)
data<-tax_glom(data,"genus")
data

data<-transform_sample_counts(data, function(x) x/sum(x))
data

# c.	Summarize the data at the genus level: number of taxa and summary and histogram of log10(mean taxa abundances)
ntaxa(data)

summary(taxa_sums(data)/nsamples(data))

hist(log10(taxa_sums(data)/nsamples(data)))

# d.	Filter out those taxa at the genus level with a mean relative abundance larger than 0.001
data<-prune_taxa((taxa_sums(data)/nsamples(data))>0.001, data)
data

# e.	Normalize the data at genus level to relative abundances
  sample_sums(data)
  data<-transform_sample_counts(data, function(x) x/sum(x))
  data

# f.	Compute alpha diversity Shannon index and the effective number of taxa as follows:
  rich<-estimate_richness(data, measures =c("Shannon"))
  effnum <- exp(rich$Shannon) 

# g.	Test for differences in the effective number of species between cases and controls (boxplot and Wilcoxon test)

  metadata<-sample_data(data)
  boxplot(rich$Shannon~metadata$disease)
  
  boxplot(effnum~metadata$disease)
  wilcox.test(rich$Shannon~metadata$disease)  

  wilcox.test(effnum~metadata$disease)  
  
# h.	Perform ordination plots using both MDS and NMDS and the BC, UF and wUF distances. Display the samples according to disease.

  BC.dist<-distance(data, method ="bray")
  UF.dist<-distance(data, method = "unifrac")
  wUF.dist<-distance(data, method = "wunifrac")
  
  BC.MDS=ordinate(data,	"MDS",	distance	=	BC.dist)
  UF.MDS=ordinate(data,	"MDS",	distance	=	UF.dist)
  wUF.MDS=ordinate(data,	"MDS",	distance	=	wUF.dist)
  
  plot_ordination(data,	BC.MDS,	color="disease")
  plot_ordination(data,	UF.MDS,	color="disease")
  plot_ordination(data,	wUF.MDS,	color="disease")

  BC.NMDS=ordinate(data,	"NMDS",	distance	=	BC.dist)
  UF.NMDS=ordinate(data,	"NMDS",	distance	=	UF.dist)
  wUF.NMDS=ordinate(data,	"NMDS",	distance	=	wUF.dist)
  
  plot_ordination(data,	BC.NMDS,	color="disease")
  plot_ordination(data,	UF.NMDS,	color="disease")
  plot_ordination(data,	wUF.NMDS,	color="disease")
  
  
# i.	Test for global differences between cases and controls with PERMANOVA

  test.adonis<-adonis2(BC.dist~metadata$disease)
  test.adonis
  
  
# 3.	Use coda4microbiome R package to analyze qin dataset at genus level
# a.	Perform an exploratory analysis of log-ratios
  
  library(coda4microbiome)
  
  x<-data.frame(otu_table(data))
  dim(x)
  x<-t(x)
  y<-metadata$disease
  lr_t2d<-explore_logratios(x,y)
  
  lr_t2d$`max log-ratio`
  lr_t2d$`names max log-ratio`
  lr_t2d$`name of most important variables`

# b.	Implement variable selection to obtain a microbial signature that is predictive of disease status

  taxanames<-colnames(x)
  colnames(x)<-paste("taxa",(1:ncol(x)), sep="")
  taxanames<-data.frame(taxanames, namestaxa=colnames(x))
  
  set.seed(123)
  
  
  t2d<-coda_glmnet(x,y)
  
  t2d$`apparent AUC`
  t2d$`mean cv-AUC`
  t2d$`sd cv-AUC`
  microbial_score<-t2d$predictions
  
  plot(metadata$BMI, microbial_score)
  
  boxplot(microbial_score~metadata$gender)
  
  plot(metadata$age, microbial_score)
  
  boxplot(microbial_score~metadata$age_category)
  
  covar=data.frame(age=metadata$age)
  
  t2d<-coda_glmnet(x,y, covar=covar)
  
  
  # set.seed(123)
   t2d<-coda_glmnet(x,y, lambda = "lambda.min")
   
  
   
  # 
  # null distribution
  
  t2dnull<-coda_glmnet_null(x,y, niter=10)  # default niter=100
  
  hist(t2dnull$accuracy,breaks = 10)
  t2dnull$`confidence interval`
  
  