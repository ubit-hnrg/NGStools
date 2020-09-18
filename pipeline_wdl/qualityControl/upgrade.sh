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
<<<<<<< HEAD
        output_path_upgrade = '/data/excecutionHNRG/upgrade_stats/'+ $2 + '/inputs'
          
        path_data_bam = '/data/resultsHNRG/'+$exp_name+'/*/*[!_haplotype].bam'
        #path_data_index = /data/resultsHNRG/'+$exp_name+'/*/*

        for f in exp_name; do readlink -f $path_data_bam; done

=======
        output_path_upgrade='/data/excecutionHNRG/upgrade_stats/'$2'/inputs'
          
        path_data='/data/resultsHNRG/'$exp_name'/*/*'
        input = $3
dirname($path_data'[!_haplotype].'$input)
>>>>>>> 899e1dfa52e11cada99b6e86e333780c742ada83


        #for f in exp_name; do readlink -f $path_data'/*[!_haplotype].bam;done

        
