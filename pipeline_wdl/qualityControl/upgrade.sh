 #!/bin/bash

 ###############33
 ## bash for upgrade_statistics; collect input mandatory files
## ##TSOxxxx
##  bam ###bam;
##  index #### bai
# usage = 
###########
        set -e

        exp_name=$1 #Nombre Experimento, Ej: TSOxxxxxx
        output_path_upgrade = '/data/excecutionHNRG/upgrade_stats/'+ $2 + '/inputs'
          
        path_data_bam = '/data/resultsHNRG/'+$exp_name+'/*/*[!_haplotype].bam'
        #path_data_index = /data/resultsHNRG/'+$exp_name+'/*/*

        for f in exp_name; do readlink -f $path_data_bam; done



        