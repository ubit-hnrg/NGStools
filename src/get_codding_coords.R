library("dplyr")
library("biomaRt")

## get all geneSymbols from agilent capture  kit SureSelect
AgilentSureselectV6_path<-"/home/ariel/repos/HNRG-pipeline-V0.1/libraries/GRCh37/S07604624_SureSelectHumanAllExonV6+UTRs_Covered_GRCh37.bed.gz"
AgilentSureselectV6 <- read.table(AgilentSureselectV6_path,header = F,sep = '\t',stringsAsFactors = F)

x = AgilentSureselectV6$V4
pat = pat = "ENST[0-9]+"
whole_ens_transcripts <-  str_extract_all(x,pat)%>%unlist%>%unique
genes <- sapply(AgilentSureselectV6$V4,function(x){strsplit(x,'\\|')[[1]][2]})%>%
    sapply(function(x){strsplit(x,'\\,')[[1]][1]})%>%unname()%>%unlist()%>%unique
    
# ensemble dataset
mart = useMart("ensembl",host="grch37.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl", mart = mart)
# print attributes and filters in the ensemble dataset
if(F){ 
    #attrib = listAttributes(ensembl)
    #filters = listFilters(ensembl)
    #write.table(filters,'~/Documents/filters_v2.tsv',sep = '\t')
    #write.table(attrib,'~/Documents/attrib_v2.tsv',sep = '\t')
}

### Utility function ####
get_canonical <- function(gene_info){
    return(gene_info[gene_info$transcript_length==max(gene_info$transcript_length),])
}                           

if(F){
    get_gene_info_by_gene <- function(genes){
        return(data.frame(unique(getBM(attributes = c("chromosome_name","genomic_coding_start","genomic_coding_end",
                                "external_gene_name","exon_chrom_start","exon_chrom_end",
                                "strand","rank","ensembl_transcript_id","transcript_length"), 
                            filters="external_gene_name", values=genes,
                            mart = ensembl))) 
        )
    }
}

get_gene_info_by_chr <- function(chr){
    return(data.frame(unique(getBM(attributes = c("chromosome_name","genomic_coding_start","genomic_coding_end",
                            "external_gene_name","exon_chrom_start","exon_chrom_end",
                            "strand","rank","ensembl_transcript_id","transcript_length","ensembl_gene_id"), 
                           filters="chromosome_name", values=chr,
                           mart = ensembl))) 
    )
}

chrs = c(1:22,'X','Y')
chrs = "5"
start_time <- Sys.time()
chr_query_subset <-  splited <- cano <- canonicos <- list()
info_by_chr <- lapply(1:length(chrs),function(i){
        query <<- get_gene_info_by_chr(chr = chrs[i]) # split in batch 
        # restrict to agilent sure select transcripts
        chr_query_subset[[i]]<- query[query$ensembl_transcript_id%in%whole_ens_transcripts,]
        splited[[i]] <<- split(chr_query_subset[[i]],f = chr_query_subset[[i]]$ensembl_gene_id)   # split individual genes
        cano[[i]] <<- lapply(splited[[i]],get_canonical)
        canonicos[[i]] <<- do.call('rbind',cano[[i]])  # batch de canonicos
        return(canonicos[[i]])
})
canoCHR <- do.call("rbind",info_by_chr)
canoCHR$external_gene_name%>%unique%>%as.character%>%length
Sys.time() - start_time
canoCHR$ensembl_gene_id%>%unique%>%as.character%>%length


genesOK <- genes[genes%in%canoCHR$external_gene_name]
revisar <- genes[!genes%in%canoCHR$external_gene_name]

ENSs <- grepl('^ENS',revisar)
LOCs <- grepl('^LOC[0-9]',revisar)
LINCs <- grepl('^LINC[0-9]',revisar)
MIRs <- grepl('^MIR',revisar)

check <- revisar[(!ENSs)&!(LOCs)&!(LINCs)&!(MIRs)]


#canonical <- read.table('/home/ariel/repos/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_canonical_transcript_GRCh37.tsv.gz',sep = '\t',header = T)
#geneSymbols <- canonical$geneSymbol%>%unique%>%as.character 
#splitedGeneSymbols<- split(geneSymbols,1:60)


#start_time <- Sys.time()
#info_by_gene <- lapply(1:length(splitedGeneSymbols),function(i){
#        gen_query_subset <- get_gene_info_by_gene(genes = splitedGeneSymbols[[i]]) # split in batch 
#        splited <- split(gen_query_subset,f =gen_query_subset$external_gene_name)   # split individual genes
#        cano <- lapply(splited,get_canonical)
#        canonicos <- do.call('rbind',cano)  # batch de canonicos
#        return(canonicos)
#})
#cano <- do.call("rbind",info_by_gene)
#cano$external_gene_name%>%unique%>%as.character%>%length
#Sys.time() - start_time







#mapp_info = data.frame(unique(getBM(attributes = c("chromosome_name",
#                            "external_gene_name","ensembl_transcript_id","transcript_length"), 
#                           #filters="ensembl_transcript_id", values="ENST00000269305",
#                           filters="external_gene_name", values="TP53",
#                           mart = ensembl)))
#mapp_info[order(mapp_info$ensembl_transcript_id),]










##################################
#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org", path="/biomart/martservice",ensemblRedirect = FALSE)
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)
#transcript_info = data.frame(unique(getBM(attributes = c("chromosome_name", "genomic_coding_start","genomic_coding_end",
#                            "external_gene_name","cds_start","cds_end","exon_chrom_start","exon_chrom_end",
#                            "strand","rank"), 
#                           filters="ensembl_transcript_id", values="ENST00000269305",
#                           mart = ensembl)))
#transcript_info[order(transcript_info$rank),]
