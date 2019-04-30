#!/usr/bin/python

import json
import pandas as pd
import os
from __future__ import division
import numpy as np
import argparse


###### importar jsons
#archivo_paths = ["/home/usuario/fastp_results/1711242.txt", "/home/usuario/fastp_results/1802140.txt","/home/usuario/fastp_results/1804642.txt","/home/usuario/fastp_results/1805817.txt","/home/usuario/fastp_results/1809341.txt","/home/usuario/fastp_results/1809720.txt","/home/usuario/fastp_results/1810255.txt","/home/usuario/fastp_results/1901364.txt","/home/usuario/fastp_results/1901981.txt","/home/usuario/fastp_results/EB761.txt","/home/usuario/fastp_results/EB790.txt","/home/usuario/fastp_results/EB802.txt"]



parser = argparse.ArgumentParser(prog='fastp_stat_report.py',description='Extract statistic info from fastp_json reports', usage='%(prog)s  --fastp_json_report --output_file')
parser.add_argument('-i','--fastp_json_report_list', help='input file with one_sample_multi_lane_fastpjson_report_file per line')
parser.add_argument('-o','--output_file', help='a tsv file with all statistics')
args = parser.parse_args()


results_dict = {}
reportes = []

json_in = args.fastp_json_report_list
out = args.output_file

 
with open(json_in) as fp:
    content = fp.readlines()
    content = [x.strip() for x in content] 

######################
for q in range(len(content)):
    array_reports = []
    total_reads_bef1 = []
    q20_rates_bef1 = [] 
    q30_rates_bef1 = [] 
    read1_mean_length_bef1 = []
    read2_mean_length_bef1 = []
    total_reads_aft = []
    q20_rates_aft = [] 
    q30_rates_aft = []
    read1_mean_length_aft = []
    read2_mean_length_aft = []
    #####filter_results
    passed_filter_reads = []
    low_quality_reads = []
    too_many_N_reads = []
    too_short_reads = []
    too_long_reads = []
    total_bases_bef1 = []
    total_bases_aft = []
    #results_dict = []
    
    lectura_media_r1_bef = 0
    lectura_media_r2_bef = 0
    lectura_media_r1_aft = 0
    lectura_media_r2_aft = 0


    
    
    lines = tuple(open(content[q], 'r'))
    lenght = len(lines)


    for i in range(lenght):

        jsonfilename=lines[i].strip('\n')
        json=pd.read_json(jsonfilename)
        
        #update observables.
        ####campos importantes before filtering
        total_reads_bef1.append(json.summary.before_filtering['total_reads'])
        total_bases_bef1.append(json.summary.before_filtering['total_bases'])
        q20_rates_bef1.append(json.summary.before_filtering['q20_rate'])
        q30_rates_bef1.append(json.summary.before_filtering['q30_rate'])
        read1_mean_length_bef1.append(json.summary.before_filtering['read1_mean_length'])
        read2_mean_length_bef1.append(json.summary.before_filtering['read2_mean_length'])
       
    
        ####campos importantes after filtering
        total_reads_aft.append(json.summary.after_filtering['total_reads'])
        total_bases_aft.append(json.summary.before_filtering['total_bases'])
        q20_rates_aft.append(json.summary.after_filtering['q20_rate'])
        q30_rates_aft.append(json.summary.after_filtering['q30_rate'])
        read1_mean_length_aft.append(json.summary.after_filtering['read1_mean_length'])
        read2_mean_length_aft.append(json.summary.after_filtering['read2_mean_length'])
      
    
        ####campos importantes filtering_results
    
        passed_filter_reads.append(json.filtering_result['passed_filter_reads'])
        low_quality_reads.append(json.filtering_result['low_quality_reads'])
        too_many_N_reads.append(json.filtering_result['too_many_N_reads'])
        too_short_reads.append(json.filtering_result['too_short_reads'])
        too_long_reads.append(json.filtering_result['too_long_reads'])
        
        sample_name=os.path.basename(content[q]).strip('.txt')
        
        
    total_lecturas=sum(total_reads_bef1)
    total_lecturas_passTrue=sum(total_reads_aft)
    porcentaje_passed = round((total_lecturas_passTrue*100)/total_lecturas,2)
    porcentaje_low_qual = round((sum(low_quality_reads)*100)/total_lecturas,2)
    porcentaje_too_many_N = round((sum(too_many_N_reads)*100)/total_lecturas,2)
    porcentaje_too_short_reads = round((sum(too_short_reads)*100)/total_lecturas,2)
    porcentaje_too_long_reads = round((sum (too_long_reads)*100)/total_lecturas,2)

        #############q20 y q30
    #porc_q30_bases_aft = round((sum(q30_bases_aft))/sum(total_bases_aft),2)
   # porc_q20_bases_aft = round((sum(q20_bases_aft)*100)/sum(total_bases_aft),2)
   # porc_q30_bases_bef1 = round((sum(q30_bases_bef1)*100)/sum(total_bases_bef),2)
   # porc_q20_bases_bef1 = round((sum(q20_bases_bef1)*100)/sum(total_bases_bef),2)

    porc_q30_rates_aft = round(np.average(q30_rates_aft,weights=total_bases_aft),2)
    porc_q20_rates_aft = round(np.average(q20_rates_aft,weights=total_bases_aft),2)
    porc_q30_rates_bef1 = round(np.average(q30_rates_bef1,weights=total_bases_bef1),2)
    porc_q20_rates_bef1 = round(np.average(q20_rates_bef1,weights=total_bases_bef1),2)
    
    #for x in range(lenght):
        #lectura_media_r1_bef += (total_reads_bef1[x]/total_lecturas)*read1_mean_length_bef1[x]
        #lectura_media_r2_bef += (total_reads_bef1[x]/total_lecturas)*read2_mean_length_bef1[x]
       # lectura_media_r1_aft += (total_reads_aft[x]/total_lecturas_passTrue)*read1_mean_length_aft[x]
        #lectura_media_r2_aft += (total_reads_aft[x]/total_lecturas_passTrue)*read2_mean_length_aft[x]
    
    lectura_media_r1_befNP = round(np.average(read1_mean_length_bef1,weights=total_bases_bef1),2)
    lectura_media_r2_befNP = round(np.average(read2_mean_length_bef1,weights=total_bases_bef1),2)
    lectura_media_r1_aftNP = round(np.average(read1_mean_length_aft,weights=total_bases_aft),2)
    lectura_media_r2_aftNP = round(np.average(read2_mean_length_aft,weights=total_bases_aft),2)
    
    
    results_dict.update({"total de lecturas antes del filtrado":int(total_lecturas)})
    results_dict.update({"total de lecturas despues del filtrado":int(total_lecturas_passTrue)})
    results_dict.update({"total de lecturas despues del filtrado[%]":porcentaje_passed})
    results_dict.update({"lecturas lowqual[%]":porcentaje_low_qual})
    results_dict.update({"lecturas NNNNN[%]":porcentaje_too_many_N})
    results_dict.update({"lecturas muy cortas[%]":porcentaje_too_short_reads})
    results_dict.update({"lecturas muy largas[%]":porcentaje_too_long_reads})
    results_dict.update({"longitud media R1 antes del filtrado":lectura_media_r1_befNP})
    results_dict.update({"longitud media R2 antes del filtrado":lectura_media_r2_befNP})
    results_dict.update({"longitud media R1 despues del filtrado":lectura_media_r1_aftNP})
    results_dict.update({"longitud media R2 despues del filtrado":lectura_media_r2_aftNP})
    results_dict.update({"q30 despues del filtrado[%]":porc_q30_rates_aft})
    results_dict.update({"q20 despues del filtrado[%]":porc_q20_rates_aft})
    results_dict.update({"q30 antes del filtrado[%]":porc_q30_rates_bef1})
    results_dict.update({"q20 antes del filtrado[%]":porc_q20_rates_bef1})

    res=pd.Series(results_dict).to_frame()
    res.columns = [sample_name]
        #res.index
    reportes.append(res.loc[[u'total de lecturas antes del filtrado',
                 u'total de lecturas despues del filtrado',
                 u'total de lecturas despues del filtrado[%]',
                 u'lecturas lowqual[%]',
                 u'lecturas NNNNN[%]',
                 u'lecturas muy cortas[%]',
                 u'lecturas muy largas[%]',
                
                 u'q30 antes del filtrado[%]',
                 u'q30 despues del filtrado[%]',
                 u'q20 antes del filtrado[%]',
                 u'q20 despues del filtrado[%]',
                             
                 u'longitud media R1 antes del filtrado',
                 u'longitud media R2 antes del filtrado',
                 u'longitud media R1 despues del filtrado',
                 u'longitud media R2 despues del filtrado'],:])
        
pd.options.display.float_format = '{:.2f}'.format
display(pd.concat(reportes, axis=1).round(2))
pd.concat(reportes, axis=1).to_csv(out,sep='\t')
