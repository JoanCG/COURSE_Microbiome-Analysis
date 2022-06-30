## PRACTICAL EXERCISE (practical microbiome.docx)

## delete existing objects
rm(list=ls())



#-------------------------------------#
#        LOAD LIBRARIES and DATA      #
#-------------------------------------#
library("phyloseq")
library("vegan")
library("coda4microbiome")

data<-readRDS("/home/jcalle/Desktop/Statistical_Methods_for_Microbiome_Analysis/qin_all.rds")



#-------------------------------------#
#               EXERCISE 2            #
#-------------------------------------#

## a)
n_taxa <- ntaxa(data)
n_samples <- nsamples(data)
n_meta <- dim(data@sam_data)[2]


## b) 
dataGenus <- tax_glom(data, "genus")
prop_dataGenus <- transform_sample_counts(dataGenus, function(x) x/sum(x))
  ## check relative abundance transformation
  sum(colSums(prop_dataGenus@otu_table))==n_samples

n_taxaGenus <- ntaxa(prop_dataGenus)
summary(rowMeans(prop_dataGenus@otu_table))
hist(log10(rowMeans(prop_dataGenus@otu_table)))


## c)
prop_dataGenus_filtered <- prune_taxa(rowMeans(prop_dataGenus@otu_table)>0.001, prop_dataGenus)
ntaxa(prop_dataGenus_filtered)
ntaxa(prop_dataGenus)


## d)
sample_sums(prop_dataGenus_filtered) # no suman 1 porque hemos filtrado los datos, HAY QUE VOLVER A NORMALIZAR
prop_dataGenus_filtered <- transform_sample_counts(prop_dataGenus, function(x) x/sum(x))
  ## check relative abundance transformation
  colSums(prop_dataGenus_filtered@otu_table)
  
  
## e)
rich<-estimate_richness(prop_dataGenus_filtered, measures =c("Shannon")) ## usamos "data" o "prop_dataGenus_filtered"?
effnum <- exp(rich$Shannon) 


## f) 
boxplot(effnum~prop_dataGenus_filtered@sam_data$study_condition)
wilcox.test(effnum~prop_dataGenus_filtered@sam_data$study_condition) # si tomamos p<0.05 como significativo, HAY DIFERENCIAS SIGNIFICATIVAS


## g) 
BC_distance <- distance(prop_dataGenus_filtered, method ="bray")
  ## MDS
  MDS_BC <- ordinate(prop_dataGenus_filtered, "MDS", distance = BC_distance)
  plot_ordination(prop_dataGenus_filtered, MDS_BC, color="disease") # NO SE VEN DIFERENCIAS
  ## NMDS
  NMDS_BC <- ordinate(prop_dataGenus_filtered, "NMDS", distance = BC_distance)
  plot_ordination(prop_dataGenus_filtered, NMDS_BC, color="disease") # NO SE VEN DIFERENCIAS

## Unifrac distance
UF_distance <- distance(prop_dataGenus_filtered, method ="unifrac")
  ## MDS
  MDS_UF <- ordinate(prop_dataGenus_filtered, "MDS", distance = UF_distance)
  plot_ordination(prop_dataGenus_filtered, MDS_UF, color="disease") # NO SE VEN DIFERENCIAS
  ## NMDS
  NMDS_UF <- ordinate(prop_dataGenus_filtered, "NMDS", distance = UF_distance)
  plot_ordination(prop_dataGenus_filtered, NMDS_UF, color="disease") # NO SE VEN DIFERENCIAS

## weighted-Unifrac distance
wUF_distance <- distance(prop_dataGenus_filtered, method ="wunifrac")
  ## MDS
  MDS_wUF <- ordinate(prop_dataGenus_filtered, "MDS", distance = wUF_distance)
  plot_ordination(prop_dataGenus_filtered, MDS_wUF, color="disease") # NO SE VEN DIFERENCIAS
  ## NMDS
  NMDS_wUF <- ordinate(prop_dataGenus_filtered, "NMDS", distance = wUF_distance)
  plot_ordination(prop_dataGenus_filtered, NMDS_wUF, color="disease") # NO SE VEN DIFERENCIAS
  
# la frontera rara que se ve en los plots (forma "triangular") es consecuencia de la composicionalidad de los datos, porque su espacio está restringido (suman 1). En la unweighted no se aprecia porque solo tiene en cuenta las relaciones filogeneticas, mientras que la WEIGHTED tiene también en cuenta las abundancias

  
## h)
permanova<-adonis2(BC_distance~prop_dataGenus_filtered@sam_data$disease)
permanova ## HAY DIFERENCIAS SIGNIFICATIVAS ENTRE LOS CASOS (T2D) Y EL CONTROL (healthy)



#-------------------------------------#
#               EXERCISE 3            #
#-------------------------------------#

## a) 
## change taxa as columns and rows as samples
prop_dataGenus_filtered@otu_table <- t(prop_dataGenus_filtered@otu_table)
  ## check if taxa are rows (result should be FALSE, "explore_logratios" function wants taxa AS COLUMNS)
  taxa_are_rows(prop_dataGenus_filtered)

## change format of data for explore_logratios function
data.x <- as.data.frame(prop_dataGenus_filtered@otu_table)
data.y <- as.factor(prop_dataGenus_filtered@sam_data$disease)

## explore logratios
logratios<-explore_logratios(x=data.x,y=data.y, measure = "glm")

logratios$`order of importance` # variables ordered by relevance
logratios$`name of most important variables` # names of relevant variables
logratios$`max log-ratio` # variable-pair with greatest logratio
logratios$`names max log-ratio` # names of variables in the pair with greatest logratio


## b)
## change colnames for plotting purposes
og_taxa <- colnames(data.x)
genus_taxa <- sapply(og_taxa, function(x) strsplit(x, "|", fixed=TRUE)[[1]][6])
colnames(data.x) <- genus_taxa

## run regression
glmnet<-coda_glmnet(x=data.x,y=data.y, lambda = exp(-2.3)) # selecciono un lambda manualmente basado en el plot, ya que el "lambda.1se" solo incluye 2 variables y no es muy informativo

sum(glmnet$`log-contrast coefficients`) # la suma de los coeficientes es 0
glmnet$taxa.name # taxas associated with the disease status
glmnet$`mean cv-AUC`

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


