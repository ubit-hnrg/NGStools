## from command line 
# conda install -c r r-essentials


###############################################################
### INSTALL OntologyIndex and   Annotation for HPO (TCGAome)###
###############################################################
install.packages(devtools)
library(devtools)

## install some dependencies
source("https://bioconductor.org/biocLite.R")
biocLite(c('AnnotationDbi''biomaRt','cowplot','GOSemSim','omicade4','ontoCAT','org.Hs.eg.db','ReactomePA','RTCGAToolbox','topGO','ellipse','mixOmics'))

install_github("ggbiplot", "vqv")
install_github('TCGAome', 'priesgo')

biocLite('ontologyIndex')
install.packages('Matrix')

###############################################################
#For jupyther kernel
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))
devtools::install_github('IRkernel/IRkernel')

###############################################################

