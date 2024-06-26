#!/usr/bin/python3
from pdf_reports import pug_to_html, write_report

import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-fq','--fastq_file', help='fasqt quality report file')
parser.add_argument('-aq','--alignment_quality', help='samtools report file')
parser.add_argument('-dq','--depth_quality', help='bams_stat_depth_global_coverage_stats')
parser.add_argument('-s','--sex', help='sex prediction file')
parser.add_argument('-n','--name', help='sample name')
parser.add_argument('-d','--date', help='date')
parser.add_argument('-t','--TSO', help='TSO name')
parser.add_argument('-o','--pdf_report') 

args =  parser.parse_args()

fastp = args.fastq_file
a_qual = args.alignment_quality
bam_dp = args.depth_quality
sex_f = args.sex
sample_id = args.name
fecha_corrida = args.date
TSO = args.TSO

outfile = args.pdf_report

###pug template
#pug_template = "/home/hnrg/HNRG-pipeline-V0.1/src/pug_template/template_ubit.pug"
pug_template = '/home/hnrg/NGStools/pipeline_wdl/unified_wes/V2.0/template_ubit.pug'

####soft hnrg_pipeline
#hnrg_soft = {'Plataforma':['Illumina Nextseq500'],'Referencia':['GRCh37.75'],'Software pipeline HNRG':['Gatk-package-4.0.8.1','fastp 0.20.0','Bwa 0.7.15-r1140','bedtools v2.28.0','SAMTOOLS_VERSION 1.9','SnpEff 4.3t','VCFtools (0.1.17)','cromwell-37'],'Anotacion':['dbSNP:All_20180423 v151','indels_site: Mills_and_1000G_gold_standard.indels.b37','hapmap_3.3','1000G_omni2.5.b37','Gwas_catalog_v1.0.2-associations_e96_r2019-04-06','ESP6500SI-V2-SSA137.snps_indels','Clinvar_20190219','PharmGKBvcf_2016','ExAC.r1.sites.vep'],'Annovar':['refGene','avsnp150','esp6500siv2_all','1000g2015aug_all','exac03','gnomad_exome','gnomad_genome','clinvar_20180603','dbscsnv11','dbnsfp35a,rmsk','tfbsConsSites','cytoBand','wgRna','targetScanS','genomicSuperDups','dgvMerged','gwasCatalog','ensGene','knownGene','intervar_20180118']}
#hnrg_pipeline = pd.DataFrame(list(hnrg_soft.values()), index=hnrg_soft.keys()).replace(np.nan, '', regex=True).T

#hnrg_pipeline = hnrg_pipeline[['Referencia','Plataforma','Software pipeline HNRG','Anotacion','Annovar']]
#hnrg_pipeline.reset_index(drop=True,inplace=True)

sex_df = pd.read_csv(sex_f, sep= '\t', header = None, usecols =[1,3,5])
sex_df.rename(columns= {1:'SEXO',3:'lecturas_X',5:'lecturas_Y'},inplace =True)

filtrado_df = pd.read_csv(fastp, sep='\t', header = None,index_col = 0).T 
antes_col = ['total de lecturas antes del filtrado','q30 antes del filtrado[%]','longitud media R1 antes del filtrado','longitud media R2 antes del filtrado']
despues_col = ['total de lecturas despues del filtrado','total de lecturas despues del filtrado[%]','q30 despues del filtrado[%]','q20 despues del filtrado[%]','longitud media R1 despues del filtrado','longitud media R2 despues del filtrado']
NNN_col = ['lecturas lowqual[%]', 'lecturas NNNNN[%]', 'lecturas muy cortas[%]','lecturas muy largas[%]']
bases_col = ['total de GIGAbases antes del filtrado','total de GIGAbases despues del filtrado','contenido GC antes','contenido GC despues']

bases = filtrado_df[bases_col]
antes = filtrado_df[antes_col]#.index('total de lecturas antes del filtrado')
desp = filtrado_df[despues_col]#.set_index('total de lecturas despues del filtrado')
desp.rename(columns = {'total de lecturas despues del filtrado':'lecturas despues','total de lecturas despues del filtrado[%]':'lecturas despues[%]','longitud media R1 despues del filtrado':'long media R1 despues','longitud media R2 despues del filtrado':'long media R2 despues'}, inplace = True)
NNN = filtrado_df[NNN_col]#.set_index('lecturas lowqual[%]')



###profundidad
prof = pd.read_csv(bam_dp, sep='\t',index_col = 0, header = None).T
col_dp = ['mean_DP','median_DP','std_DP','q25_DP','q75_DP','q95_DP']
dp_col = ['dp>=1','dp>=10','dp>=20','dp>=30','dp>=50','dp>=100']
stats_dp = prof[col_dp]
dp = prof[dp_col]

###alineamiento
alineamiento = pd.read_csv(a_qual,sep='\t',index_col =0,header = None).T
reads_col = ['Number-of-reads-properly-paired','Percent-of-reads-properly-paired-in-Library','reads-duplicated-in-Library','reads-MQ0-in-Library']
rates_col = ['error-rate-in-Library','efficienciy-in-Library','filtering_ratio','on_target_ratio-in-Library','unduplicated_ratio-in-Library']
Gbases_col = ['GBases_before_trimm','GBases_after_trimm','Gbases-mapped-(cigar)-in-Library','Gbases-without-dup-mapped-in-Library']
lib_col = ['Target-size-(GBPs)','maximum-length-in-Library']

al_reads = alineamiento[reads_col]
al_rates = alineamiento[rates_col]
al_Gbases = alineamiento[Gbases_col]
al_lib = alineamiento[lib_col]


html = pug_to_html(pug_template,title=sample_id, genero=sex_df, TSO = TSO, fecha = fecha_corrida, bases = bases, antes = antes, desp = desp, NNN = NNN,dp = dp, stats_dp = stats_dp, al_reads = al_reads, al_rates = al_rates, al_Gbases = al_Gbases, al_lib = al_lib)#, profundidad= muestra_prof_lib,version= version, N_lecturas_antes = N_lecturas_antes, N_lecturas_dps = N_lecturas_dps, N_lecturas_dps_percent = N_lecturas_dps_percent, lecturas_lowqual = lecturas_lowqual, lecturas_NN = lecturas_NN, short_reading = short_reading, too_long=too_long, q_30_bef=q_30_bef,q_30_aft=q_30_aft, q20_bef= q20_bef, q20_after=q20_after)
write_report(html, outfile)

