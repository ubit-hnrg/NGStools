 #!/bin/bash
        set -e

        exp_name = $1 #Nombre Experimento, Ej: TSOxxxxxx
        output_path_upgrade = '/data/excecutionHNRG/upgrade_stats/'+ $2 + '/inputs'
          
        path_data = '/data/resultsHNRG/'+$exp_name+'/*/*'

        for f in exp_name; do readlink -f $path_data'/*[!_haplotype].bam;done

        