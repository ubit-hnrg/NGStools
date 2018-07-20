library('ontologyIndex')
data(hpo)
header = hpo = read.csv("/home/ariel/HPO/archive/annotation/ALL_SOURCES_FREQUENT_FEATURES_phenotype_to_genes.txt",nrows = 1,sep='/t',header = F)

header = c('')
hpo = read.csv("/home/ariel/HPO/archive/annotation/ALL_SOURCES_FREQUENT_FEATURES_phenotype_to_genes.txt",nrows = 10,sep='\t',skip = 1,header = F)
