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
        #output_path_upgrade='/data/excecutionHNRG/upgrade_stats/'$2'/inputs'
          
        path_bam='/data/resultsHNRG/'$exp_name'/*/*[!_haplotype].bam'
        #input = $3
        #dirname($path_data'[!_haplotype].'$input)


        for f in exp_name; do readlink -f $path_bam;done

        
