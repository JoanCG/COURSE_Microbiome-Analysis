library(coda4microbiome)
library(mia)
library(dplyr)
library(ggpubr)
library(phyloseq)
# Get the data:
data<-readRDS("hansen.rds")
data

### Explore phyloseq data
ntaxa(data)
nsamples(data)
sample_sums(data)

# Focus on the metadata: number of subjects and number of time points for each subject.

metadata<-data@sam_data

unique(metadata$subject_id)
length(unique(metadata$subject_id))

    # days_from_first_collection variable exploring

table(sample_data(data)$subject_id)

table(sample_data(data)$subject_id, metadata$Diet)


hist(sample_data(data)$days_from_first_collection)

ggscatter(data=data.frame(sample_data(data)), 
          x="days_from_first_collection",
          y= "subject_id",
        color="subject_id"
        )

ggscatter(data=data.frame(sample_data(data)), 
          y="days_from_first_collection",
          x= "subject_id",
          color="subject_id"
)

    # Agglomerate at Species level
psq_s <- tax_glom(data, "species") 


    # Summary data at species level
ntaxa(psq_s)
nsamples(psq_s)
sample_sums(psq_s)

    # Normalize
sp_prop <- transform_sample_counts(psq_s, function(x) {x/sum(x)})
sample_sums(sp_prop)

summary(taxa_sums(sp_prop)/nsamples(sp_prop))
hist(taxa_sums(sp_prop)/nsamples(sp_prop))

    # Filter >0.0001
sp_filt <- filter_taxa(sp_prop, function (x) mean(x)> 0.0001, TRUE)
sp_filt

    # Normalize filtered data
sp_filt_prop <- transform_sample_counts(sp_filt, function(x) {x/sum(x)})
sp_filt_prop

### Apply coda4microbiome longitudinal
otutab <- data.frame(otu_table(sp_filt))  # Extract abundance matrix
metadata <- data.frame(sample_data(sp_filt)) # Extract metadata

# Order data
o<-order(metadata$ID_num)

metadata<-metadata[o,]

### Set function parameters:
x <- impute_zeros(t(otutab)) # microbiome abundance

x<-x[o,]

x_time <- metadata$days_from_first_collection  # observation times
subject_id <- metadata$ID_num  # subject id
y <- as.factor(metadata$Diet)  # diet 
ini_time = 0
end_time = 280

## Apply coda_glmnet_longitudinal():
lr.hansen <- explore_lr_longitudinal(x,y, x_time, subject_id, ini_time, end_time)
lr.hansen

dftaxanames<-data.frame(taxa=paste("taxa",(1:ncol(x)),sep=""),taxaname=colnames(x))

colnames(x)<-paste("taxa",(1:ncol(x)),sep="")

## Apply coda_glmnet_longitudinal():
set.seed(123)

hansen <-coda_glmnet_longitudinal (x,y, x_time, subject_id, ini_time, end_time, lambda="lambda.min",nfolds=4)
  # Results:
hansen$`predictions plot`
hansen$`signature plot`
hansen$`trajectories plot`
hansen$`apparent AUC`
hansen$`mean cv-AUC`
hansen$taxa.name
  