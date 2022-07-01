## PRACTICAL EXERCISE (practical microbiome.docx)

## delete existing objects
rm(list=ls())



#-------------------------------------#
#        LOAD LIBRARIES and DATA      #
#-------------------------------------#
library("phyloseq")
library("vegan")
library("coda4microbiome")

data<-readRDS("/home/jcalle/Desktop/Statistical_Methods_for_Microbiome_Analysis/day4/hansen.rds")



#-------------------------------------#
#               EXERCISE 2            #
#-------------------------------------#

## a)
ntaxa(data)
nsamples(data)
dim(data@sam_data)[2]


## b) 
length(unique(data@sam_data$subject_id))
length(unique(data@sam_data$days_from_first_collection))


## c)
hist(data@sam_data$days_from_first_collection)
range(data@sam_data$days_from_first_collection)

## d)
data <- tax_glom(data, "species")
data <- transform_sample_counts(data, function(x) x/sum(x))
## check relative abundance transformation
colSums(data@otu_table)
  
  
## e)
ntaxa(data)
hist(log10(rowMeans(data@otu_table)))

## f) 


## g) 
data_filtered <- prune_taxa(rowMeans(data@otu_table)>0.0001, data)
ntaxa(data_filtered)
ntaxa(data)
  
## h)
data_filtered <- transform_sample_counts(data_filtered, function(x) x/sum(x))
## check relative abundance transformation
colSums(data_filtered@otu_table)



#-------------------------------------#
#               EXERCISE 3            #
#-------------------------------------#

## a) 
## change taxa as columns and rows as samples
data_filtered@otu_table <- t(data_filtered@otu_table)
  ## check if taxa are rows (result should be FALSE, "explore_logratios" function wants taxa AS COLUMNS)
  taxa_are_rows(data_filtered)
  




## change format of data and define parameters for explore_logratios function
x <- as.data.frame(data_filtered@otu_table) # microbiome abundances
metadata <- data_filtered@sam_data
x_time = metadata$days_from_first_collection;  # observation times
subject_id = metadata$ID_num; # subject id
y= as.factor(metadata$Diet); # diet
ini_time = 0;
end_time = 280;

## order metadata and OTUs according to ID_num
metadata <- metadata[order(metadata$ID_num),]
x <- x[rownames(metadata),]

## explore logratios_longitudidanl
logratios<-explore_lr_longitudinal(x,y, x_time, subject_id, ini_time, end_time)

logratios$`order of importance` # variables ordered by relevance
logratios$`name of most important variables` # names of relevant variables
logratios$`max log-ratio` # variable-pair with greatest logratio
logratios$`names max log-ratio` # names of variables in the pair with greatest logratio


## b)
## set a seed to replicate results
set.seed(123)

## change colnames for plotting purposes
og_taxa <- colnames(x)
species_taxa <- sapply(og_taxa, function(x) strsplit(x, "|", fixed=TRUE)[[1]][7])
colnames(x) <- species_taxa

ecam <-coda_glmnet_longitudinal(x,y, x_time, subject_id, ini_time, end_time, lambda="lambda.min",nfolds=4)

sum(ecam$`log-contrast coefficients`) # la suma de los coeficientes es 0
ecam$taxa.name # taxas associated with the disease status
ecam$`apparent AUC`
ecam$`mean cv-AUC`


## test significance of our prediction with null distribution
glmnet_null <- coda_glmnet_null(data.x,data.y, niter = 10)

glmnet_null$`confidence interval`

hist(glmnet_null$accuracy, breaks = 10)


# # las distribuciones tan extrañas en la clasificación (segundo plot) se pueden deber a que la variable de interés ("disease") esté correlacionada con otra variable. HABRÍA QUE CORREGIR EL MODELO
# boxplot(glmnet$predictions ~ prop_dataGenus_filtered@sam_data$age) # parece que aquí puede haber algo
# 
#   ## remove NAs to run corrected glmnet
#   data_filtered.x <- data.x[-which(is.na(prop_dataGenus_filtered@sam_data$age)),]
#   data_filtered.y <- data.y[-which(is.na(prop_dataGenus_filtered@sam_data$age))]
#   covar <- prop_dataGenus_filtered@sam_data$age[-which(is.na(prop_dataGenus_filtered@sam_data$age))]
#   
# glmnet_corrected <- coda_glmnet(x=data_filtered.x, y=data_filtered.y, covar = covar)


