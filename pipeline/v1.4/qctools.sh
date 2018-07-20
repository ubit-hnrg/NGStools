#!/bin/bash

## Script version
VERSION="1.4"

## DECLARE VARIABLES


declare -A tools=(["BWA"]="/home/bitgenia/tools/aligners/bwa-0.7.12" ["BT2"]="/home/bitgenia/tools/aligners/bowtie2-2.2.6" ["BBMAP"]="/home/bitgenia/tools/aligners/bbmap" ["PICARD"]="/home/bitgenia/tools/accessoryTools/picard-tools-2.0.1"  ["GATK"]="/home/bitgenia/tools/variantCallers" ["PRINSEQ"]="/home/bitgenia/tools/accessoryTools/prinseq-lite-0.20.4" ["PLATYPUS"]="/home/bitgenia/tools/variantCallers/Platypus_0.8.1" ["SNPEFF"]="/home/bitgenia/tools/functionalAnn/snpEff" ["BPLATANN"]="/home/bitgenia/tools/accessoryTools" ["FASTQC"]="/home/bitgenia/tools/accessoryTools/FastQC");


############### FUNCTIONS #######################

function trim {
    trimmed=$1
    trimmed=${trimmed%% }
    trimmed=${trimmed## }

    echo "$trimmed"
}



function createResultDirectories {
    LOG=$RESULTPATH/$SAMPLENAME\_logs
    RESULT=$RESULTPATH/$SAMPLENAME\_Result
    TIMELOG=$RESULTPATH/$SAMPLENAME\_time_logs
    PREFIX_QC=$SAMPLENAME\_QC

    mkdir -p $LOG   # make directory for log files
    mkdir -p $TIMELOG
    mkdir -p $RESULT
}

function getDataToTrim {
    echo -n "Please enter the number of position or length to trim:"
    read POSLEN
    getInteractiveSampleInfo
}

function getDataToGCRange {
    echo -n "Please enter the GC content range (50-60,75-90):"
    read RANGE
    getInteractiveSampleInfo
}

function getDataComplete {
   echo -n "Please enter the number of position or length to trim:"
   read POSLEN
   echo -n "Please enter the GC content range (50-60,75-90):"
   read RANGE
   getInteractiveSampleInfo
}


function getInteractiveSampleInfo {

    echo -n "Please enter the sample name:"
    read SAMPLENAME

    
    SAMPLENAME=$(trim ${SAMPLENAME}) 
    if [ -z "$SAMPLENAME" ]
     then
       echo "Error: Sample name is null or empty"
       exit
     fi

    SM="$SAMPLENAME"  # RGSM=String    REQUIRED sample name Required.
    PREFIX="$SAMPLENAME"
    echo -n "Please enter fastq file containing the 1st reads from the pair (absolute path /home/bitgenia/...):"
    read FIRSTREADS

   
    FIRSTREADS=$(trim ${FIRSTREADS}) 
    if [ -z "$FIRSTREADS" ] || [ ! -f "$FIRSTREADS" ]
     then
        echo "Error: fastq file is null or not exists"
        exit
     fi

    echo -n "Please enter fastq file containing mate pairs (absolute path /home/bitgenia/...):"
    read MATES
    
    MATES=$(trim ${MATES})
    if [ -z "$MATES" ] || [ ! -f "$MATES" ]
     then
        echo "Error: fastq file is null or not exists"
        exit
     fi

    echo -n "Please enter results directory (absolute path /home/bitgenia/...):"
    read RESULTPATH
   
    RESULTPATH=$(trim ${RESULTPATH})
    if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
     then
        echo "Error: the result directory is null or not exists"
        exit
     fi

    createResultDirectories

}

function fastqc {
   echo "$(date) Running Sample $SM, FastQC"
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv ${tools["FASTQC"]}/fastqc -t 6 $FIRSTREADS $MATES --outdir=$RESULT  > $LOG/$PREFIX\_fastqc.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: FastQC command"
      exit
    fi
}

function trimming {
    
      if [ -n "${1}" ]
       then
            CMD=${1};
      fi

    echo "$(date) Running Sample $SM, PrinSEQ-lite"
    #Quality control. Input: 2 fastq files (one per mate); output: 2 fastq files corrected with qc methods
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  perl ${tools["PRINSEQ"]}/prinseq-lite.pl -fastq $FIRSTREADS -fastq2 $MATES -out_good $RESULT/$PREFIX\_TRIMMED -out_bad $RESULT/$PREFIX\_TRIMMED_BAD ${CMD} ${POSLEN} -log $LOG/$PREFIX\_prinseq.log > $LOG/$PREFIX\_prinseq_results.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: PrinSEQ command"
      exit
    fi
}
function gc_content {
    
    echo "$(date) Running Sample $SM, PrinSEQ-lite"
    #Quality control. Input: 2 fastq files (one per mate); output: 2 fastq files corrected with qc methods
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  perl ${tools["PRINSEQ"]}/prinseq-lite.pl -fastq $FIRSTREADS -fastq2 $MATES -out_good $RESULT/$PREFIX\_GC -out_bad $RESULT/$PREFIX\_GC_BAD -range_gc ${RANGE} -log $LOG/$PREFIX\_prinseq.log > $LOG/$PREFIX\_prinseq_results.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: PrinSEQ command"
      exit
    fi
}

function completeTrimmingGC {
    
      if [ -n "${1}" ]
       then
            CMD=${1};
      fi

    echo "$(date) Running Sample $SM, PrinSEQ-lite"
    #Quality control. Input: 2 fastq files (one per mate); output: 2 fastq files corrected with qc methods
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  perl ${tools["PRINSEQ"]}/prinseq-lite.pl -fastq $FIRSTREADS -fastq2 $MATES -out_good $RESULT/$PREFIX\_FIX -out_bad $RESULT/$PREFIX\_FIX_BAD ${CMD} ${POSLEN} -range_gc ${RANGE} -log $LOG/$PREFIX\_prinseq.log > $LOG/$PREFIX\_prinseq_results.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: PrinSEQ command"
      exit
    fi
}

function interactiveScript {

    while(true)
    do
        echo
        printf "Choose from the following operations: \n"
        printf "[F]astQC\n"
        printf "[T]rimming\n"
        printf "[G]C_content\n"
        printf "[A]All\n"
        read -p "Your choice: " choice
        echo   "-----------------------------------------"
        case $choice in
        [fF])
          getInteractiveSampleInfo
          fastqc
          ;;
        [tT])
            printf "Choose from trimming command:\n"
            printf "[a]trim_to_len\n"
            printf "[b]trim_left\n"
            printf "[c]trim_right\n"
            read -p "Your choice: " choicePS
          
            case $choicePS in
              [a])
                getDataToTrim
                trimming "-trim_to_len"    
              ;;
              [b])
                getDataToTrim
                trimming "-trim_left"
              ;;
              [c])
                getDataToTrim
                trimming "-trim_right"
              ;;
              *)
                 echo "wrong choice!"
                 exit;
            esac
        ;;
        [gG])
             getDataToGCRange
             gc_content
            ;;
        [aA])
            printf "Choose from trimming command:\n"
            printf "[a]trim_to_len\n"
            printf "[b]trim_left\n"
            printf "[c]trim_right\n"
            read -p "Your choice: " choicePS
          
            case $choicePS in
              [a])
                getDataComplete
                completeTrimmingGC "-trim_to_len"    
              ;;
              [b])
                getDataComplete
                completeTrimmingGC "-trim_left"
              ;;
              [c])
                getDataComplete
                completeTrimmingGC "-trim_right"
              ;;
              *)
                 echo "wrong choice!"
                 exit;
            esac
            ;;
        *)
           echo "wrong choice!"
           exit;
        esac
        read -p "Do you wish to continue? [y]es or [n]o: " ans
        ans=$(trim ${ans})
        if [ $ans == 'n' ]
            then
             echo "Exiting the script. Have a nice day!"
             break
            else 
             continue
        fi
    done
}


########################### MAIN ######################
 echo
 echo   "---------------------------------------------"
 echo   "   Welcome to QC tools ${VERSION}"
 echo   "   For support and documentation contact: info@bitgenia.com"
 echo   "   www.bitgenia.com - www.biargentina.com.ar"
 echo   "   Copyright (c) 2013-2016 BITGENIA - BIA"
 echo   "---------------------------------------------"

interactiveScript
