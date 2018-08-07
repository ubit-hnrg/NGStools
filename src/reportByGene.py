#!/usr/bin/python

import argparse
import os
import pip
from pyensembl import EnsemblRelease as ensbl
import numpy as np
import pandas as pd
import mygene


import io
import time
import subprocess
import sys
def popen_with_logging(cmd,logfile = 'out.log'):
    with io.open(logfile, 'wb') as writer, io.open(logfile, 'rb', 1) as reader:
        process = subprocess.Popen(cmd, stdout=writer)
        while process.poll() is None:
            sys.stdout.write(reader.read())
            time.sleep(0.5)
        # Read the remaining
        sys.stdout.write(reader.read())
    return None


def write_and_logging(mje,writer=None,stdout = True):
    if writer is None:
        writer = io.open('out.log', 'wb')
    if stdout:
        print mje
    writer.write(mje)

# ## Utility functions

def parseArgs(description = 'Get coverage at different DPs for a given list of genes.'):
    parser = argparse.ArgumentParser(description= description)
    parser.add_argument('-b','--bamfile',required=True,help='input bam file')
    parser.add_argument('-g','--genelistfile',required=True)
    parser.add_argument('-o','--outpath',required=True)
    parser.add_argument('-l','--logfile',required=False,default= 'out.log')
    parser.add_argument('-r','--ref',dest = 'ref',default = '',required = True,type=int,help = 'Reference genome is mandatory since only a list of genames is given and we need to match it inside our bam file')
    parser.add_argument('-v','--vcffile',default = '',required=False,help='input vcf file to be usde for complementary analysis', )
    parser.add_argument('-p','--prefix',dest = 'prefix',default = '',required = False)
    parser.add_argument('-s','--splitBamFile', dest='split', action='store_false',help = 'in next iterations this flag will allow to split the reducedBAM file by chromosomes')
    parser.add_argument('-e','--entrez', dest='entrez', action='store_true',help = 'use entrez ids')
    parser.set_defaults(split=False)
    parser.set_defaults(entrez=False)

    args = parser.parse_args()
    # MANDATORY ARGUMENTS
    bamfile = args.bamfile
    vcffile = args.vcffile
    outpath = args.outpath
    prefix = args.prefix
    ref = args.ref
    genelistfile = args.genelistfile
    logfile = args.logfile    
    #OPTIONAL ARGUMENTS
    split= args.split
    entrez= args.entrez
    
    return bamfile, vcffile, outpath, prefix, ref, genelistfile, split, entrez, logfile
    
    
    
def check_ref(ref,outpath):
    if ref not in [19,37,38]:
        print 'please, reference must be one of 19,37,38\n'
        print 'bye'
        quit()
    else:
        if ref == 19:
            print 'using ref hg19, '
        if ref == 37:
            print 'using ref hgrc37'
        if ref == 38:
            print 'using ref hgrc38'
    print ''
    
    if not os.path.isdir(outpath):
        os.system('mkdir %s'%outpath)    

    
def get_genelist(genelistfile,entrez,writer = None):
    genes = pd.read_csv(genelistfile,header = None)
    if not entrez:
        genelist = genes.iloc[:,0].tolist()
        mje =  'read genes from file:\n'
        write_and_logging(mje,writer = writer)
        write_and_logging(' '.join(genelist),writer = writer)
        #print '\n'
    if entrez:
        genelist = genes.iloc[:,1].tolist()
        mje =  'read genes from file:'
        write_and_logging(mje,writer = writer)
        write_and_logging(' '.join(genelist),writer = writer)
        #print '\n'
    return genelist
    
    
def test_mode():
    ref = 37
    genelistfile = './test/genelist'
    bamfile = './test/test_chr2_190-192k.bam'
    outpath = './test_output/'
    split = False
    return(ref,genelistfile,bamfile,outpath,split)
    
    
def get_locis(genelist,ref,outpath,entrez = False,write_bedfile = True,writer =None):
    analyzed_genes = []
    try:
        if (ref == 37)|(ref == 19):
            ensemble_version = 75
    except:
        'You must to install the ensemble version 75, please run \n'        
        'pyensembl install --release 75 --species homo_sapiens \n'
    try:
        if ref == 38:
            ensemble_version = 87
    except:
        'You must to install latest ensemble version (87), please run \n'        
        'pyensembl install --release 87 --species homo_sapiens \n'
        
    data = ensbl(ensemble_version)
    locis = []
    mg = mygene.MyGeneInfo()
    
    if entrez == False:    
        for gen in genelist:
            gene = gen.strip(' ').replace('-', '')
            try:
                loc =  data.loci_of_gene_names(gene)
                locis.append(loc[0].to_dict())
                analyzed_genes.append(gene)
            except ValueError:
                gene_query = mg.query(gene, size=5,species = 'human')
                gids = gene_query['hits']
                if(len(gids)==0):
                    mje = 'warning, gene %s was not found and ignored\n'%gene 
                    write_and_logging(mje,writer = writer)
                    continue
                else:
                    j=0
                    for i in range(len(gids)): 
                        gid = int(gids[i][u'_id'])
                        try:
                            #print gid
                            ensblID = mg.getgene(gid, fields='ensembl',species = 'human')[u'ensembl'] [u'gene']
                            trylocus = data.locus_of_gene_id(ensblID)
                            break
                        except:
                            j=j+1
                            continue

                    # check if there was a match 
                    if i == j:
                        mje = 'warning, gene %s was not found and ignored\n'%gene 
                        write_and_logging(mje,writer = writer)
                        continue
                        
                    if len(gids)>1:
                        mje = 'warning: more than one mathch with %s'%gene
                        write_and_logging(mje,writer = writer)
                        mje = 'gen id: %s (ensmbl id %s ) was assumed'%(gid,ensblID)
                        write_and_logging(mje,writer = writer)
                loc =  data.locus_of_gene_id(ensblID)
#                print gene 
#                print loc
#                locis.append(loc[0].to_dict())
                locis.append(loc.to_dict())
                analyzed_genes.append(gene)

    if entrez == True:
        for gene in genelist:
            try:
                ensblID = mg.getgene(gene, fields='ensembl',species = 'human')[u'ensembl'] [u'gene']
                loc =  data.locus_of_gene_id(ensblID)
                locis.append(loc[0].to_dict())
                analyzed_genes.append(gene)
            except:
                mje = 'warning, gene %s was not found and ignored\n'%gene 
                write_and_logging(mje,writer = writer)
                continue
                
    
    
    print '\n'
    ggenes = analyzed_genes
    locis = pd.DataFrame(locis)
    locis['name'] = ggenes
    locis.rename(columns = {'contig':'chrom'},inplace = True)
    locis['score'] = 0 ## this is only for completness, in order to recognize strand as the right next field
    locis = locis[['chrom','start','end','name','score','strand']]
    if ref == 19:
        locis['chrom'] = locis.chrom.apply(lambda x: 'chr'+str(x))

    # write file
    bedfile = outpath+'gene_loci.bed'
    if write_bedfile:
        locis.to_csv(bedfile,sep = '\t',index = False,header = None)
        
    return locis, bedfile, locis.columns


def run_samtools_view(bamfile,bedfile,outpath,split=False):
    base = os.path.basename(bamfile).split('.bam')[0]
    outbam = outpath + base + '_reduced.bam'
    if not split:
        bedfile = bedfile
        call1 = 'samtools view -bh -L %s %s > %s'%(bedfile,bamfile,outbam)
    status = os.system(call1)    
    if status ==0:
        print 'ok samtools'
    else:
        print 'samtools have failed'
        #quit()

    # get index
    call2 = 'samtools index %s'%(outbam)
    index = os.system(call2)
    if index ==0:
        print 'ok index'
    else:
        print 'index have failed'
        #quit()

    return  outbam
    

def serie_counting(x):
    return pd.Series({'bp_at_10dp':(x>10).sum(),'bp_at_20dp':(x>20).sum(),'bp_at_30dp':(x>30).sum()})

def run_bedtools_coverage(bam,bed,outpath,exonbed = False,prefix =''):
    if exonbed :
        resolution = 'exon'
    else:
        resolution = 'gene'

    outcoveragefile = outpath +'.'+prefix + resolution +'.tsv' # lo dejo como archivo oculto
    absolute_file_report =  outpath + '.'+ prefix + resolution +'_absolute_report'+'.tsv'  # lo dejo como archivo oculto
    relative_file_report= outpath + prefix + resolution +'_relative_report'+'.tsv'

    print outcoveragefile
    call = 'coverageBed -a %s -b %s > %s'%(bed,bam,outcoveragefile)
    status = os.system(call)
    if status ==0:
        print 'ok coveragedBed'
    else:
        print 'coveragedBed have failed'

    results = pd.read_table(outcoveragefile,header = None)
    results.columns = ['chr','start','end','name','score','strand','matched_reads','coverage_dp1','len_gen_bp','rel_coverage_dp1']
        
    
    ### the later approach give us the coverage as the 
    #number of read bases at least ONE time over the total length of our gene of interest
    
    ## in the following I atempt to compute the coverage at different minimal dp levels, namely: at dp= 10, 20, 30
    
    #compute dp per base x gene
    file_by_position = outpath + resolution +'_coverage_by_position.txt'
    call2 = 'coverageBed -a %s -b %s -d > %s'%(bed,bam,file_by_position)
    os.system(call2)
    
    bypos = pd.read_csv(file_by_position,sep = '\t',header = None)

    grouping = bypos.groupby([3])[7].apply(lambda x: serie_counting(x))  
    
    res  = grouping.unstack(); 
    #reset index and rename it
    cols = res.columns  #reserve columns
    res.reset_index(inplace = True)
    res.columns = ['name']+list(cols)

    
    ###  merge both analysis
    coverage = pd.merge(results,res, on = 'name')
    
    ### relativise results
    abscols = ['bp_at_10dp','bp_at_20dp','bp_at_30dp']



    
    if not exonbed:
        coverage.to_csv(absolute_file_report,index = False, sep = '\t')
        relatives = np.round((coverage[abscols].transpose() /coverage['len_gen_bp'].values).transpose(),7)
        coverage[abscols] = relatives
        coverage.to_csv(relative_file_report,index = False, sep = '\t')


    if exonbed:
        coverage.rename(columns = {'score':'Gene'}).to_csv(absolute_file_report,index = False, sep = '\t')
        relatives = np.round((coverage[abscols].transpose() /coverage['len_gen_bp'].values).transpose(),7)
        coverage[abscols] = relatives
        coverage.rename(columns = {'score':'Gene'}).to_csv(relative_file_report,index = False, sep = '\t')


    if False:    
#        exome_report = coverage.groupby(['name'])[['coverage_dp1','len_gen_bp','bp_at_10dp','bp_at_20dp','bp_at_30dp']].sum()

        exome_report = coverage.groupby(['name'])[['len_gen_bp','bp_at_10dp','bp_at_20dp','bp_at_30dp']].sum()
        for col in abscols:
            exome_report[col] =  exome_report[col]/exome_report['len_gen_bp']
        exome_coverage = exome_report.copy()
        #rel = exome_report.apply(lambda x:x[abscols]/float(x['len_gen_bp']),axis = 1)#.drop([u'len_gen_bp'],axis = 1)
        #print rel.head()
        #exon_coverage = pd.merge(exome_report[['len_gen_bp']]),rel,how = 'inner',left_index = True, right_index = True)
        exome_coverage.reset_index(inplace = True)
        #exome_report.reset_index(inplace = True)
        exome_coverage.rename(columns = {'score':'Gene'},inplace = True)

        abs_coverage = coverage.copy()
        abs_coverage[exome_report.columns] = exome_report
        abs_coverage.to_csv(absolute_file_report,index =False , sep = '\t')

        rel_coverage = coverage.copy()
        rel_coverage[exome_coverage.columns] = exome_coverage
        rel_coverage.to_csv(relative_file_report,index =False , sep = '\t')

    
    
    return outcoveragefile


def get_exons(genelist,ref,outpath,entrez = False, write_bedfile = True,writer = None):
    analyzed_genes = []
    try:
        if (ref == 37)|(ref == 19):
            ensemble_version = 75
    except:
        'You must to install the ensemble version 75, please run \n'        
        'pyensembl install --release 75 --species homo_sapiens \n'
    try:
        if ref == 38:
            ensemble_version = 87
    except:
        'You must to install latest ensemble version (87), please run \n'        
        'pyensembl install --release 87 --species homo_sapiens \n'
        
    data = ensbl(ensemble_version)
    ############ By Exons #############
    exonbed = []
    mg = mygene.MyGeneInfo()

    if entrez == False:
        for gen in genelist:
            gene = gen.strip(' ').replace('-', '')
            try:
                exons = data.exon_ids_of_gene_name(gene)
                exonLocus = [data.locus_of_exon_id(e) for e in exons]
                exonLoci = [ex.to_dict() for ex in exonLocus]
                exonLoci = pd.DataFrame(exonLoci)
                exonLoci['name'] = exons
                exonLoci['score'] = [gene]*len(exons)
                exonLoci.rename(columns = {'contig':'chrom'},inplace = True)
                #exonLoci['score'] = 0 ## this is only for completness, in order to recognize strand as the right next field
                exonLoci = exonLoci[['chrom','start','end','name','score','strand']]
                exonbed.append(exonLoci)
                analyzed_genes.append(gene)
            except ValueError:
                gene_query = mg.query(gene, size=5,species = 'human')
                gids = gene_query['hits']
                if(len(gids)==0):
                    mje = 'warning, gene %s was not found and ignored\n'%gene 
                    write_and_logging(mje,writer = writer)
                    continue
                else:
                    j=0
                    for i in range(len(gids)): 
                        gid = int(gids[i][u'_id'])
                        try:
                            ensblID = mg.getgene(gid, fields='ensembl',species = 'human')[u'ensembl'] [u'gene']
                            trylocus = data.locus_of_gene_id(ensblID)
                            break
                        except ValueError:
                            j=j+1
                            continue

                    # check if there was a match 
                    if i == j:
                        mje = 'warning, gene %s was not found and ignored\n'%gene 
                        write_and_logging(mje,writer = writer)
                        continue
                        
                    if len(gids)>1:
                        mje =  'warning: more than one mathch with %s'%gene
                        write_and_logging(mje,writer = writer)
                        mje = 'gen id: %s (ensmbl id %s ) was assumed'%(gid,ensblID)
                        write_and_logging(mje,writer = writer)
                exons = data.exon_ids_of_gene_id(ensblID)
                exonLocus = [data.locus_of_exon_id(e) for e in exons]
                exonLoci = [ex.to_dict() for ex in exonLocus]
                exonLoci = pd.DataFrame(exonLoci)
                exonLoci['name'] = exons
                exonLoci['score'] = [gene]*len(exons)
                exonLoci.rename(columns = {'contig':'chrom'},inplace = True)
                #exonLoci['score'] = 0 ## this is only for completness, in order to recognize strand as the right next field
                exonLoci = exonLoci[['chrom','start','end','name','score','strand']]
                exonbed.append(exonLoci)
                analyzed_genes.append(gene)

    if entrez == True:
        for gene in genelist:
            try:
                ensblID = mg.getgene(gene, fields='ensembl',species = 'human')[u'ensembl'] [u'gene']                
                exons = data.exon_ids_of_gene_id(ensblID)
                exonLocus = [data.locus_of_exon_id(e) for e in exons]
                exonLoci = [ex.to_dict() for ex in exonLocus]
                exonLoci = pd.DataFrame(exonLoci)
                exonLoci['name'] = exons
                exonLoci['score'] = [gene]*len(exons)
                exonLoci.rename(columns = {'contig':'chrom'},inplace = True)
                #exonLoci['score'] = 0 ## this is only for completness, in order to recognize strand as the right next field
                exonLoci = exonLoci[['chrom','start','end','name','score','strand']]
                exonbed.append(exonLoci)
                analyzed_genes.append(gene)
            except ValueError:
                mje = 'warning, gene %s was not found and ignored\n'%gene 
                write_and_logging(mje,writer = writer)
                continue
                

    print '\n'
    exonbed = pd.concat(exonbed)
    if ref == 19:
        exonbed.chrom = exonbed.chrom.apply(lambda x: 'chr'+str(x))
    
    # write file
    exonbedfile = outpath+'exons_loci.bed'

    if write_bedfile:
        exonbed.to_csv(exonbedfile,sep = '\t',index = False,header = None)
    
        
    return exonbed, exonbedfile, exonbed.columns

def main(test = False):
    if test:
        ref, genelistfile, bamfile, outpath, split = test_mode()
    else:
        bamfile, vcffile, outpath, prefix, ref, genelistfile, split, entrez, logfile = parseArgs()
    #check params
    check_ref(ref=ref,outpath = outpath)

    wlogfile = io.open(logfile, 'wb')

    genelist = get_genelist(genelistfile,entrez=entrez,writer = wlogfile)
    gene_loci , bedfile , locicolumns = get_locis(genelist,ref = ref,outpath=outpath,entrez=entrez,writer = wlogfile)
    exon_loci , exon_bedfile , exon_columns = get_exons(genelist,ref = ref,outpath=outpath,entrez = entrez,writer = wlogfile)

    
    reduced_bamfile = run_samtools_view(bamfile,bedfile,split=split,outpath = outpath)
    
    #compute coverage by gen along the bamfile
    coverage_file = run_bedtools_coverage(reduced_bamfile,bedfile,outpath,prefix =prefix)
    coverage_by_exon = run_bedtools_coverage(reduced_bamfile,exon_bedfile,outpath,exonbed= True,prefix = prefix)
    
    
    #os.system('rm %s'%exon_bedfile)
    wlogfile.close()
    #os.system('rm %s'%gene_loci)
    #os.system('rm %s'%exon_loci)



if __name__ == "__main__":
    main()

