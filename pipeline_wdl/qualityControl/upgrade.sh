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
        exp_lib=$2
        trim_front=$3
        trim_tail=$4

        #output_path_upgrade='/data/excecutionHNRG/upgrade_stats/'$2'/inputs'
          
        path_bam='/data/resultsHNRG/'$exp_name'/*/*[!_haplotype].bam'
        path_index='/data/resultsHNRG/'$exp_name'/*/*[!_haplotype].bai'
        path_excecution='/home/hnrg/executionsHNRG/upgrade_stats/'$exp_name
        path_fasq='/data/excecutionsHNRG/'$exp_name'/inputs/'$exp_name'.txt'
        #input = $3
        #dirname($path_data'[!_haplotype].'$input)

        mkdir -p $path_excecution'/upgrade_stats/inputs'
        for f in exp_name; do readlink -f $path_index >> $path_excecution/inputs/bams_index.txt; readlink -f $path_bam >> $path_excecution/inputs/bams.txt; cp $path_fasq $path_excecution/inputs; done


        python json_upgrade.py -sa $path_excecution'/inputs/'$1'.txt' -bam $path_excecution'/inputs/bams.txt' -bai $path_excecution'/inputs/bams_index.txt' -tf $3  -tt $4 -e $2 
