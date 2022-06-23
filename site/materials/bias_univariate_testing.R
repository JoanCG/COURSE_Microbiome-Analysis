##########################################################
#
# BIAS UNIVARIATE TESTING WITH LOG-PROPORTIONS AND CLR
# Toy example
#
#########################################################


library(ggpubr)
library(tidyr)
library(MASS)

set.seed(1)

# Binary outcome: 

y<-c(rep(0,100),rep(1,100))
y<-factor(y, labels=c("controls","cases"))

# Example 1: 5 taxa equally abundant (0.2,0.2,0.2,0.2,0.2)
# counts for each sample are randomly generated according to a multinomial distribution 

M<-  rmultinom(200, 500, c(0.2,0.2,0.2,0.2,0.2))

# Example 2: 5 taxa with mean abundances (0.02,0.02,0.02,0.47,0.47)
# counts for each sample are randomly generated according to a multinomial distribution 

# M<-  rmultinom(200, 500, c(0.02,0.02,0.02,0.47,0.47))

# Transpose abundance table so that rows are samples and columns are taxa

M<-t(M)

# Fold-change for each taxa for cases with respect to controls:

F<-c(5,1/2,1/10,1,1)

# Multiply the abundances of cases by F: 

M[(101:200),]<-t(apply(M[(101:200),],1,function(x) x*F))

# Normalization to obtain relative abundances:

x<-M/rowSums(M)
rowSums(x)

# Data-frame with log(relative abundances) and outcome:

dflogx<-data.frame(log(x),y)
colnames(dflogx)<-c("logX1","logX2","logX3","logX4","logX5","y")

# clr-transformation:

clrx<-log(x) # CLR se calcula sobre la MATRIZ DE ABUNDANCIAS RELATIVAS
clrx<-clrx-(rowSums(log(x))/ncol(x))
rowSums(clrx) # como clr centra los datos, la suma de las abundancias debe ser 0


# Data-frame with clr-abundances and outcome:

dfclrx<-data.frame(clrx,y)
colnames(dfclrx)<-c("clrX1","clrX2","clrX3","clrX4","clrX5","y")

# Box-plot log(relative abundances) vs outcome:

ggboxplot(dflogx, x = "y", y = c("logX1","logX2","logX3","logX4","logX5"),
          color = "y", palette = "npg",
          #combine = TRUE,
          merge = "flip",
          add = "point", shape = "Grade", ylab=FALSE)

# Wilcoxon test for log(relative abundances):

wilcox.test(dflogx[,1]~y)
wilcox.test(dflogx[,2]~y)
wilcox.test(dflogx[,3]~y)
wilcox.test(dflogx[,4]~y)
wilcox.test(dflogx[,5]~y)

# Box-plot clr-abundances vs outcome:

ggboxplot(dfclrx, x = "y", y = c("clrX1","clrX2","clrX3","clrX4","clrX5"),
          color = "y", palette = "npg",
          #combine = TRUE,
          merge = "flip",
          add = "point", shape = "Grade", ylab=FALSE)


# Wilcoxon test for clr-abundances:

wilcox.test(dfclrx[,1]~y)
wilcox.test(dfclrx[,2]~y)
wilcox.test(dfclrx[,3]~y)
wilcox.test(dfclrx[,4]~y)
wilcox.test(dfclrx[,5]~y)

