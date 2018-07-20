#!/bin/bash

## Script version
VERSION="1"
UPDATE="27/02/2018"

## DECLARE VARIABLES
USR='sample'
PASS='sample321654'
SERVER=192.168.210.77

declare -A tools=(["BWA"]="/home/bitgenia/tools/aligners/bwa-0.7.16a" ["BT2"]="/home/bitgenia/tools/aligners/bowtie2-2.2.6" ["BBMAP"]="/home/bitgenia/tools/aligners/bbmap" ["PICARD"]="/home/bitgenia/tools/accessoryTools/picard-2.11" ["GATK"]="/home/bitgenia/tools/variantCallers" ["PRINSEQ"]="/home/bitgenia/tools/accessoryTools/prinseq-lite-0.20.4" ["PLATYPUS"]="/home/bitgenia/tools/variantCallers/Platypus_0.8.1" ["SNPEFF"]="/home/bitgenia/tools/functionalAnn/snpEff" ["BPLATANN"]="/home/bitgenia/tools/accessoryTools" ["FASTQC"]="/home/bitgenia/tools/accessoryTools/FastQC");


declare -A GRCh37=(["HGREF"]="/home/bitgenia/bundle/hgref" ["BT2_REF"]="/home/bitgenia/bundle/hgref/human_g1k_v37_decoy" ["REF"]="/home/bitgenia/bundle/hgref/human_g1k_v37_decoy.fasta" ["KGINDELS"]="/home/bitgenia/bundle/hgref/1000G_phase1.indels.b37.vcf" ["MILLSINDELS"]="/home/bitgenia/bundle/hgref/Mills_and_1000G_gold_standard.indels.b37.vcf" ["DBSNP"]="/home/bitgenia/bundle/hgref/dbsnp_150.vcf.gz" ["GWASCAT"]="/home/bitgenia/bundle/hgref/gwascatalog.txt"  ["GWASCAT_ADJUSTED"]="/home/bitgenia/bundle/hgref/gwascatalog_adjusted.vcf" ["EVS"]="/home/bitgenia/bundle/hgref/ESP6500SI-V2-SSA137.snps_indels.vcf" ["CADD"]="/home/bitgenia/bundle/hgref/CADD_1000G.vcf" ["CLNVAR"]="/home/bitgenia/bundle/hgref/clinvar_20170801_fixed.vcf.gz" ["PHARMGKB"]="/home/bitgenia/bundle/hgref/PharmGKBvcf_2016.vcf" ["ExAC"]="/home/bitgenia/bundle/hgref/ExAC.r0.3.1.sites.vep.vcf.gz" ["HAPMAP"]="/home/bitgenia/bundle/hgref/hapmap_3.3.b37.vcf" ["KOMNI"]="/home/bitgenia/bundle/hgref/1000G_omni2.5.b37.vcf" ["KNOWNVARS"]="/home/bitgenia/scripts/known_variants_analysis-devel/grch37mybaby525.vcf" ["DBNSFP"]="/home/bitgenia/bundle/hgref/dbNSFP2.9.txt.gz" ["GRCH_VERSION"]="GRCh37.75" ["REGIONS_MT"]="/home/bitgenia/bundle/hgref/library/mitochondrial.bed");
# ojo, cambie a 2.9 temporalmente hasta tener la 2.9.2

declare -A hg19=(["HGREF"]="/home/bitgenia/bundle/hg19" ["BT2_REF"]="/home/bitgenia/bundle/hg19" ["REF"]="/home/bitgenia/bundle/hg19/ucsc.hg19.fasta" ["KGINDELS"]="/home/bitgenia/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf" ["MILLSINDELS"]="/home/bitgenia/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" ["DBSNP"]="/home/bitgenia/bundle/hg19/dbsnp_144.hg19.vcf.gz" ["GWASCAT"]="/home/bitgenia/bundle/hg19/gwascatalog.txt" ["GWASCAT_ADJUSTED"]="/home/bitgenia/bundle/hg19/gwascatalog_adjusted.vcf" ["EVS"]="/home/bitgenia/bundle/hg19/ESP6500SI-V2-SSA137.snps_indels.vcf" ["CADD"]="/home/bitgenia/bundle/hg19/CADD_1000G.vcf" ["CLNVAR"]="/home/bitgenia/bundle/hg19/clinvar_20170801_fixed_hg19.vcf.gz" ["PHARMGKB"]="/home/bitgenia/bundle/hg19/PharmGKBvcfhg19_2016.vcf" ["ExAC"]="/home/bitgenia/bundle/hg19/ExAC.r0.3.1.sites.vep.vcf.gz" ["HAPMAP"]="/home/bitgenia/bundle/hg19/hapmap_3.3.hg19.sites.vcf" ["KOMNI"]="/home/bitgenia/bundle/hg19/1000G_omni2.5.hg19.sites.vcf" ["KNOWNVARS"]="/home/bitgenia/scripts/known_variants_analysis-devel/clinvar_Result/hg19mybaby525.vcf" ["DBNSFP"]="/home/bitgenia/bundle/hg19/dbNSFP2.9.txt.gz" ["GRCH_VERSION"]="hg19" ["REGIONS_MT"]="/home/bitgenia/bundle/hg19/library/mitochondrial.bed");

#XMX=Java heap max size
#NT= Num of data threads sent to processor, Http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
#NCT= Num of CPU threads for each data thread, http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
#QFILTER_READS PARAMETER OPTIONAL min.qual.mean: "Minimum mean quality" TYPE INTEGER(Filter reads with quality score mean below the given value.)
#MINPRUN minPruning sets the kmer size for the BRUIJN tree in variant calling.
#CALL Variant Call Confidence when creating .vcf file.
#EMIT  Variant Emition Confidence when creating .vcf file.
declare -A params=(["XMX"]="16g" ["NT"]="11" ["NCT"]="8" ["QFILTER_READS"]="10" ["MINPRUN"]="3" ["CALL"]="30.0" ["EMIT"]="10.0");


#LB= RGLB=String   REQUIRED Library.
#PL=     # RG PL=String  REQUIRED platform (e.g. illumina, solid). IMPORTANTE para GATK: no puede ponerse UNKNOWN por mÃ¡s que lo muestre como posible. Debe especificarse una plataforma, la marca en general (ILLUMINA y no HISEQ200, por ejemplo).
#PU=     # RGPU=String   REQUIRED platform unit (eg. run barcode). {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
#ID=1    # RGID=String   ID Default value: 1. This option can be set to 'null' to clear the default value.
#CN=null # RGCN=String   sequencing center name Default value: null.
#DS=null # RGDS=String   description Default value: null.
#DT=null # RGDT=Iso8601Date      run date Default value: null.
#PI=null # RGPI=INTEGER  predicted insert size Default value: null.
declare -A readGroups=(["LB"]="1" ["PL"]="ILLUMINA" ["PU"]="unknown" ["ID"]="1" ["CN"]="" ["DS"]="" ["DT"]="" ["PI"]="");

#[design ID]_Padded.bed - This BED file contains a single track of the genomic regions that you can expect to sequence when using the design for target enrichment. To determine these regions, the program extends the regions in the Covered BED file by 100 bp on each side.


declare -A libraries=();

declare -A references=();

declare -A optional_params=();

############### FUNCTIONS #######################

function trim {
    trimmed=$1
    trimmed=${trimmed%% }
    trimmed=${trimmed## }

    echo "$trimmed"
}

function getReferences {

    if [ -n "${1}" ] && [ -n "${2}" ]
       then
         choice=${1};
         choiceK=${2};
       else
            printf "Choose from the following genome reference:\n"
            printf "[a] GRch37\n"
            printf "[b] hg19\n"
            read -p "Your choice: " choice
            echo   "-----------------------------------------"
      fi
      echo

      declare -a arr=()
      case $choice in
          [aA])
            for K in "${!GRCh37[@]}"; do references[$K]=${GRCh37[$K]}; done
            arr=(/home/bitgenia/bundle/hgref/library/*.bed);
          ;;
          [bB])
            for K in "${!hg19[@]}"; do references[$K]=${hg19[$K]}; done
            arr=(/home/bitgenia/bundle/hg19/library/*.bed)
          ;;
          *)
             echo "wrong choice!"
             exit;
      esac

      if [ -z $choiceK ]
        then
          printf "Choose from the following capture prep Kit:\n"
          idx=0;
          for bedpath in ${arr[@]}; do
            filename=${bedpath##*/}
            printf "[$idx] $filename\n"
            ((idx++))
          done
          read -p "Your choice: " choiceK
          echo   "-----------------------------------------"
      fi

      tLen=${#arr[@]}
      if [ ${choiceK} -lt ${tLen} ]
        then
            filename=$(basename "${arr[$choiceK]}")
            path=$(dirname "${arr[$choiceK]}")
            extension="${filename##*.}"
            filename="${filename%.*}"
            libraries["REGIONS"]=${arr[$choiceK]}
            libraries["REGIONSTXT"]=$path"/auxiliary/"$filename"_platypus.txt"
            libraries["REGIONSFB"]=$path"/auxiliary/"$filename"_freebayes.bed"
            libraries["REGIONSCOVERAGE"]=$path"/auxiliary/"$filename"_coverage.bed"
         else
            echo "wrong bedfile selected!"
            exit;
      fi
      printf 'We will use the bed file %s\n' "${libraries['REGIONS']}"
}

function createResultDirectories {
    LOG=$RESULTPATH/$SAMPLENAME\_logs
    RESULT=$RESULTPATH/$SAMPLENAME\_Result
    TIMELOG=$RESULTPATH/$SAMPLENAME\_time_logs

    mkdir -p $LOG   # make directory for log files
    mkdir -p $TIMELOG
    mkdir -p $RESULT
}

function getNonInteractiveSampleInfo {

    SAMPLENAME=$(trim ${1})
    if [ -z "${SAMPLENAME}" ]
     then
       echo "Error: Sample name is null or empty"
       exit
     fi

    SM="${SAMPLENAME}"  # RGSM=String    REQUIRED sample name Required.
    PREFIX="${SAMPLENAME}"
    FIRSTREADS=$(trim ${2})
    if [ -z "$FIRSTREADS" ] || [ ! -f "$FIRSTREADS" ]
     then
        echo "Error: fastq file is null or not exists"
        exit
     fi

    MATES=$(trim ${3})
    if [ -z "$MATES" ] || [ ! -f "$MATES" ]
     then
        echo "Error: fastq file is null or not exists"
        exit
     fi

    RESULTPATH=$(trim ${4})
    if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
     then
        echo "Error: the result directory is null or not exists"
        exit
     fi

    createResultDirectories

    for i in $(trim ${5}) $(trim ${6}); do
      KEY=${i%=*};
      VAL=${i#*=};
      optional_params[$KEY]=$VAL
    done

    DOMAIN=${optional_params["domain"]}
    if [[ -z "${DOMAIN}" ]]
     then
        echo "Warning: without domain id we can't not upload vcf file to B_Platform's ftp"
        DOMAIN=0
    fi

    BEDFILE=${optional_params["bed"]}
    if [[ -z "${BEDFILE}" ]]
     then
         BEDFILE=1 #By default use /home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV5_Padded.bed
    fi
    getReferences "a" ${BEDFILE}
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

    SM="$SAMPLENAME"
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

function processLanes {
  #FIRSTREADS es un txt con un fastq por linea
  #MATES es un txt con un fastq por linea, del 2do elemento del par en el mismo orden que firstreads
  mapfile -t R1_FILES < $FIRSTREADS
  mapfile -t R2_FILES < $MATES
  
  for (( i=0; i < ${#R1_FILES[@]}; i++))
    do
      LANE_NUMBER=$((i+1))
      FIRSTREADS=${R1_FILES[$i]}
      MATES=${R2_FILES[$i]}
      preprocessing
    done
  
  preprocessing_join
}

function preprocessing {

    PREFIX_LANE=$SAMPLENAME\_L$LANE_NUMBER
    PREFIX_LANE_QC=$SAMPLENAME\_L$LANE_NUMBER\_QC

    echo "$(date) Running Sample $SM, FastQC - Lane $LANE_NUMBER"
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX_LANE\_time.log.csv ${tools["FASTQC"]}/fastqc -t 6 $FIRSTREADS $MATES --outdir=$RESULT  > $LOG/$PREFIX_LANE\_fastqc.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: FastQC command"
      exit
    fi
    
    echo "$(date) Running Sample $SM, PrinSEQ-lite - Lane $LANE_NUMBER"
    #Quality control. Input: 2 fastq files (one per mate); output: 2 fastq files corrected with qc methods
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX_LANE\_time.log.csv  perl ${tools["PRINSEQ"]}/prinseq-lite.pl -fastq $FIRSTREADS -fastq2 $MATES -out_good $RESULT/$PREFIX_LANE_QC -out_bad $RESULT/$PREFIX_LANE_QC\_BAD -min_qual_mean ${params["QFILTER_READS"]} -log $LOG/$PREFIX_LANE\_prinseq.log > $LOG/$PREFIX_LANE\_prinseq_results.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: PrinSEQ command"
      exit
    fi

    #2018-03-19 - Arreglar los ReadGroups https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
    echo "$(date) Running Sample $SM, BWA: mem - Lane $LANE_NUMBER"
    OUTPUTFILE=$RESULT/$PREFIX_LANE\_aligned.sam
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX_LANE\_time.log.csv ${tools["BWA"]}/bwa mem -M -R "@RG\tID:${readGroups[ID]}\tSM:$SM\tLB:${readGroups[LB]}\tPL:${readGroups[PL]}" -t ${params["NT"]} ${references["REF"]} $RESULT/$PREFIX_LANE_QC\_1.fastq $RESULT/$PREFIX_LANE_QC\_2.fastq > ${OUTPUTFILE} 2> $LOG/$PREFIX_LANE\_bwa_mem.log

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: BWA command"
      exit
    fi

    echo "$(date) Running Sample $SM, PICARD: CleanSam - Lane $LANE_NUMBER"
    INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX_LANE\_cleaned.sam
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX_LANE\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar CleanSam I=${INPUTFILE} O=${OUTPUTFILE} CREATE_INDEX=true > $LOG/$PREFIX_LANE.cleansam.log 2>&1
    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard CleanSam"
      exit
    fi

    echo "$(date) Running Sample $SM, PICARD: SortSam - Lane $LANE_NUMBER"
    INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX_LANE\_sorted.bam
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX_LANE\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar SortSam I=${INPUTFILE} O=${OUTPUTFILE} SO=coordinate CREATE_INDEX=true > $LOG/$PREFIX_LANE.sortsam.log 2>&1

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard SortSam"
      exit
    fi

    echo "$(date) Running Sample $SM, SamTools: FlagStat - Lane $LANE_NUMBER"
    samtools flagstat ${OUTPUTFILE} > $RESULT/$PREFIX_LANE\_PreMarkDuplicatesStats.txt 2>&1
    
    ## Crear lista de bams con extension ".list"
    echo ${OUTPUTFILE} >> $RESULT/$SM\_bam.list
}

function preprocessing_join {

##---------- 2018-03-19 - Hasta aca per lane

    #https://gatkforums.broadinstitute.org/gatk/discussion/9588/dynamically-pass-multiple-input-to-picards-markduplicates-multiplexed-data
    awk '{print "INPUT="$1}' $RESULT/$SM\_bam.list > $RESULT/$SM\_bam.tmp
    readarray -t bam_array < $RESULT/$SM\_bam.tmp

    ## MarkDuplicates unifica los lanes https://software.broadinstitute.org/gatk/documentation/article.php?id=3060
    echo "$(date) Running Sample $SM, PICARD: MarkDuplicates"
    INPUTFILE=$(echo ${bam_array[@]})
    OUTPUTFILE=$RESULT/$PREFIX\_marked.bam
    #Marcar duplicados. Input: bam, output: bam.
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false AS=true MAX_RECORDS_IN_RAM=10000000 ${INPUTFILE} O=${OUTPUTFILE} M=$RESULT/$PREFIX\_dedupped.metrics.txt CREATE_INDEX=true > $LOG/$PREFIX.markduplicates1.log 2>&1

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard markduplicates"
      exit
    fi

##---------- 2018-03-19 - Desde aca sigue todo igual

    echo "$(date) Running Sample $SM, SamTools: FlagStat"
    samtools flagstat ${OUTPUTFILE} > $RESULT/$PREFIX\_PostMarkDuplicatesStats.txt 2>&1

    echo "$(date) Running Sample $SM, PICARD: CollectAlignmentSummaryMetrics"
    INPUTFILE=${OUTPUTFILE}
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar CollectAlignmentSummaryMetrics I=${INPUTFILE} O=$RESULT/$PREFIX\_collectAlignmentSummaryMetrics.txt R=${references["REF"]} CREATE_INDEX=true > $LOG/$PREFIX.collectAlignmentSummaryMetrics.log 2>&1

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard CollectAlignmentSummaryMetrics"
      exit
    fi


    #BQSR from the Broadâs GATK allows you to reduce the effects of analysis artefacts produced by your sequencing machines. It does this in two steps, the first analyses your data to detect covariates and the second compensates for those covariates by adjusting quality scores.

    echo "$(date) Running Sample $SM, GATK: BaseRecalibrator"
    INPUTFILE=${OUTPUTFILE}
    #https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["GATK"]}/GenomeAnalysisTK.jar -T BaseRecalibrator -I ${INPUTFILE} -R ${references["REF"]} -knownSites ${references["KGINDELS"]} -knownSites ${references["MILLSINDELS"]} -knownSites ${references["DBSNP"]} -nct ${params["NCT"]} -o $RESULT/$PREFIX.recal.table >$LOG/$PREFIX.baserecalibrator.log 2>&1


    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: GATK: BaseRecalibrator"
      exit
    fi

    echo "$(date) Running Sample $SM, GATK: PrintReads"
    INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX\.bam
    #https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T PrintReads -R ${references["REF"]} -I ${INPUTFILE} -BQSR $RESULT/$PREFIX.recal.table -nct ${params["NCT"]}  -o ${OUTPUTFILE} > $LOG/$PREFIX.printreads.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: GATK: PrintReads"
      exit
    fi
    INPUTFILE=${OUTPUTFILE}
    
    preprocessingCoverage&

}

function preprocessingCoverage {
   echo -e "#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean" > $RESULT/$PREFIX\_depthAnalysis.csv
    for region in $(cat ${libraries["REGIONSCOVERAGE"]} | sed 's/\t/:/' | sed 's/\t/-/'); do
      echo -n "$region" | tr :- \\t
      samtools depth "${INPUTFILE}" -a -r $region | awk 'BEGIN {min=100000000; max=0;start=0;s=0};{c++;s+=$3;if(start==0)start=$2;if($3>max)max=$3;if($3<=min)min=$3;l=$2-start}END {print "\011",((l>0)?l:0),"\011",((min<100000000)?min:0),"\011",max,"\011",s/((c>0)?c:1)}'
    done  >> $RESULT/$PREFIX\_depthAnalysis.csv;
    echo -e "#chrom\tstart\tend\t#readsQueCubrenA\t#basesCubiertasDeA\tsizeOfA\t%OfACubiertas" > $RESULT/$PREFIX\_breadthAnalysis.csv
    bedtools coverage -a ${libraries["REGIONSCOVERAGE"]} -b ${INPUTFILE} >> $RESULT/$PREFIX\_breadthAnalysis.csv;
}

function getInteractiveBamInfo {
   if [ ! -f "$INPUTFILE" ]
     then
        echo -n "Please enter the sample name:"
        read SAMPLENAME

        SAMPLENAME=$(trim ${SAMPLENAME})
        PREFIX="$SAMPLENAME"
        SM="$SAMPLENAME"
        if [ -z "$SAMPLENAME" ]
        then
          echo "Error: Sample name is null or empty"
          exit
        fi
        echo -n "Please enter bam file (absolute path from home):"
        read INPUTFILE

        INPUTFILE=$(trim ${INPUTFILE})
        if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
        then
          echo "Error: bam file is null or not exists"
          exit
        fi
        echo -n "Please enter results directory (full path /home/bitgenia/...):"
        read RESULTPATH

        RESULTPATH=$(trim ${RESULTPATH})
        if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
         then
            echo "Error: the result directory is null or not exists"
            exit
        fi
        createResultDirectories

        if [ -n "$1" ]
        then
            echo -n "Please enter vcf file with known variants (absolute path from home):"
            read KNOWNVARS

            KNOWNVARS=$(trim ${KNOWNVARS})
            if [ -z "$KNOWNVARS" ] || [ ! -f "$KNOWNVARS" ]
            then
                echo "Error: variants file is null or not exists"
                exit
             else
                references["KNOWNVARS"]=${KNOWNVARS};
            fi
        fi
    fi
}

function getInteractiveVcfInfo {
   if [ ! -f "$INPUTFILE" ]
     then
        echo -n "Please enter the sample name:"
        read SAMPLENAME

        SAMPLENAME=$(trim ${SAMPLENAME})
        PREFIX="$SAMPLENAME"
        if [ -z "$SAMPLENAME" ]
        then
          echo "Error: Sample name is null or empty"
          exit
        fi
        echo -n "Please enter vcf file (full path from home):"
        read INPUTFILE
        arr=( "${INPUTFILE}")
        for filepath in ${arr[@]}; do
          filepath=$(trim ${filepath})
          if [ -z "$filepath" ] || [ ! -f "$filepath" ]
          then
            echo "Error: $filepath file is null or not exists"
            exit
          fi
        done

        echo -n "Please enter results directory (full path /home/bitgenia/...):"
        read RESULTPATH

        RESULTPATH=$(trim ${RESULTPATH})
        if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
         then
            echo "Error: the result directory is null or not exists"
            exit
        fi

        createResultDirectories
    fi
}

function variantDiscovery {
    if [ -n "${1}" ]
       then
            choice=${1};
       else
            printf "Choose from the following Variant Caller:\n"
            printf "[s]amtools\n"
            printf "[f]reebayes\n"
            printf "[h]aplotypeCaller\n"
            printf "[p]latypus\n"
            read -p "Your choice: " choice
            echo   "-----------------------------------------"
    fi

    OUTPUTFILE=$RESULT/$PREFIX\_variant_calling.vcf
    case $choice in
        [sS])
            echo "$(date) Running Sample $SM, Samtools: mpileup & bcftools call"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv samtools mpileup -d8000 -uf ${references["REF"]} -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR ${INPUTFILE} -l ${libraries["REGIONS"]} | bcftools call  -mv -Ov  > ${OUTPUTFILE} 2> $LOG/$PREFIX\_samtools.log
            if [ "$?" -ne 0 ]
                then
                    echo "## ERROR: samtools mpileup & bcftools call"
                    exit
            fi
         ;;
        [fF])
            echo "$(date) Running Sample $SM, freebayes"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv freebayes  -t ${libraries["REGIONSFB"]} -f ${references["REF"]} -b ${INPUTFILE} | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & DP > 20" -F "LowQual"> ${OUTPUTFILE} 2> $LOG/$PREFIX\_freebayes.log
            if [ "$?" -ne 0 ]
                then
                    echo "## ERROR: freebayes"
                    exit
            fi
         ;;
        [hH])
            echo "$(date) Running Sample $SM, GATK: HaplotypeCaller"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T HaplotypeCaller -L ${libraries["REGIONS"]} -L ${references["REGIONS_MT"]} -nct ${params["NCT"]} -I ${INPUTFILE} --dbsnp ${references["DBSNP"]} -o $RESULT/$PREFIX\_haplotypecaller.vcf -stand_call_conf ${params["CALL"]} -stand_emit_conf ${params["EMIT"]} --minPruning ${params["MINPRUN"]} > $LOG/$PREFIX.haplotypecaller.log 2>&1
            #You are correct that 1 exome is not enough data for VQSR. Instead of VQSR, you can try hard flitering. Please read about hard filtering here: http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set
            #If you want to use VQSR, you will first need to get more data from the 1000Genomes data webpage here: http://www.1000genomes.org/data
            #You should get at least 30 or more bam files from samples that are genetically similar to your sample exome and run all of the steps involved in doing a joint analysis. Please read about how to do a joint analysis here: http://www.broadinstitute.org/gatk/guide/article?id=3893

            echo "$(date) Running Sample $SM, GATK: Hard Filtering"
            #https://www.broadinstitute.org/gatk/guide/article?id=3225
            #http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

            #echo "$(date) Running Sample $SM, GATK: SelectVariants snps"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T SelectVariants -R ${references["REF"]} -nt  ${params["NT"]} -V $RESULT/$PREFIX\_haplotypecaller.vcf -selectType SNP -o $RESULT/$PREFIX\_raw_snps.vcf > $LOG/$PREFIX\_selectVariant.raw_snps.log 2>&1

            #undefined variable ReadPosRankSum undefined variable MQRankSum
            #echo "$(date) Running Sample $SM, GATK: VariantFiltration snps"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T VariantFiltration -R ${references["REF"]} -V $RESULT/$PREFIX\_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || MQ0 >= 2 || QUAL < 30 || DP < 10 || GQ < 30" --filterName "LowQual" -o $RESULT/$PREFIX\_filtered_snps.vcf > $LOG/$PREFIX\_VariantFiltration.raw_snps.log 2>&1

            #echo "$(date) Running Sample $SM, GATK: SelectVariants indels"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T SelectVariants -R ${references["REF"]} -nt ${params["NT"]} -V $RESULT/$PREFIX\_haplotypecaller.vcf -selectType INDEL -o $RESULT/$PREFIX\_raw_indels.vcf > $LOG/$PREFIX\_selectVariant.raw_indels.log 2>&1

            #echo "$(date) Running Sample $SM, GATK: VariantFiltration indels"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T VariantFiltration -R ${references["REF"]} -V $RESULT/$PREFIX\_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "LowQualindels" -o $RESULT/$PREFIX\_filtered_indels.vcf  > $LOG/$PREFIX\_VariantFiltration.raw_indels.log 2>&1

            #echo "$(date) Running Sample $SM, GATK: CombineVariants"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["GATK"]}/GenomeAnalysisTK.jar -T CombineVariants -R ${references["REF"]} -nt ${params["NT"]} --assumeIdenticalSamples --variant $RESULT/$PREFIX\_filtered_snps.vcf --variant $RESULT/$PREFIX\_filtered_indels.vcf -o $RESULT/$PREFIX\_concatenated.vcf -genotypeMergeOptions UNIQUIFY > $LOG/$PREFIX\_combinevariants.log 2>&1

            #echo "$(date) Running Sample $SM, sed"
            sed "s/\t${SM}.variant2//" $RESULT/$PREFIX\_concatenated.vcf > $RESULT/$PREFIX\_removesamplevariant2.vcf 2> $LOG/$PREFIX\_final1_sed.log

            sed "s/\t${SM}.variant/\t${SM}/" $RESULT/$PREFIX\_removesamplevariant2.vcf > ${OUTPUTFILE} 2> $LOG/$PREFIX\_final2_sed.log
         ;;
        [pP])
            echo "$(date) Running Sample $SM, platypus"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv python ${tools["PLATYPUS"]}/Platypus.py callVariants --nCPU ${params["NT"]} --bamFiles=${INPUTFILE} --refFile=${references["REF"]}  --regions=${libraries["REGIONSTXT"]} --output=${OUTPUTFILE} > $LOG/$PREFIX\_platypus.log 2>&1
           if [ "$?" -ne 0 ]
            then
                echo "## ERROR: platypus"
                exit
           fi
         ;;
        *)
           echo "wrong choice!"
           exit;
    esac
    INPUTFILE=${OUTPUTFILE}
}

function fullvariantDiscovery {

            OUTPUTFILE=$RESULT/$PREFIX\_variant_calling.vcf
            #HaplotypeCaller
            echo "$(date) Running Sample $SM, GATK: HaplotypeCaller"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T HaplotypeCaller -L ${libraries["REGIONS"]} -L ${references["REGIONS_MT"]} -nct ${params["NCT"]} -I ${INPUTFILE} --dbsnp ${references["DBSNP"]} -o $RESULT/$PREFIX\_haplotypecaller.vcf -stand_call_conf ${params["CALL"]} -stand_emit_conf ${params["EMIT"]} --minPruning ${params["MINPRUN"]} > $LOG/$PREFIX.haplotypecaller.log 2>&1
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T SelectVariants -R ${references["REF"]} -nt  ${params["NT"]} -V $RESULT/$PREFIX\_haplotypecaller.vcf -selectType SNP -o $RESULT/$PREFIX\_raw_snps.vcf > $LOG/$PREFIX\_selectVariant.raw_snps.log 2>&1

            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T VariantFiltration -R ${references["REF"]} -V $RESULT/$PREFIX\_raw_snps.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || MQ0 >= 2 || QUAL < 30 || DP < 10 || GQ < 30" --filterName "LowQual" -o $RESULT/$PREFIX\_filtered_snps.vcf > $LOG/$PREFIX\_VariantFiltration.raw_snps.log 2>&1

            echo "$(date) Running Sample $SM, GATK: SelectVariants indels"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T SelectVariants -R ${references["REF"]} -nt ${params["NT"]} -V $RESULT/$PREFIX\_haplotypecaller.vcf -selectType INDEL -o $RESULT/$PREFIX\_raw_indels.vcf > $LOG/$PREFIX\_selectVariant.raw_indels.log 2>&1

            echo "$(date) Running Sample $SM, GATK: VariantFiltration indels"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T VariantFiltration -R ${references["REF"]} -V $RESULT/$PREFIX\_raw_indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "LowQualindels" -o $RESULT/$PREFIX\_filtered_indels.vcf  > $LOG/$PREFIX\_VariantFiltration.raw_indels.log 2>&1

            echo "$(date) Running Sample $SM, GATK: CombineVariants"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["GATK"]}/GenomeAnalysisTK.jar -T CombineVariants -R ${references["REF"]} -nt ${params["NT"]} --assumeIdenticalSamples --variant $RESULT/$PREFIX\_filtered_snps.vcf --variant $RESULT/$PREFIX\_filtered_indels.vcf -o $RESULT/$PREFIX\_concatenated.vcf -genotypeMergeOptions UNIQUIFY > $LOG/$PREFIX\_combinevariants.log 2>&1

            echo "$(date) Running Sample $SM, sed"
            sed "s/\t${SM}.variant2//" $RESULT/$PREFIX\_concatenated.vcf > $RESULT/$PREFIX\_removesamplevariant2.vcf 2> $LOG/$PREFIX\_final1_sed.log

            sed "s/\t${SM}.variant/\t${SM}/" $RESULT/$PREFIX\_removesamplevariant2.vcf > $RESULT/$PREFIX\_haplotypefinal.vcf 2> $LOG/$PREFIX\_final2_sed.log
         
            ####platypus
            echo "$(date) Running Sample $SM, platypus"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv python ${tools["PLATYPUS"]}/Platypus.py callVariants --nCPU ${params["NT"]} --minReads 10 --bamFiles=${INPUTFILE} --refFile=${references["REF"]}  --regions=${libraries["REGIONSTXT"]} --output=$RESULT/$PREFIX\_platypus.vcf --filterDuplicates=0 > $LOG/$PREFIX\_platypus.log 2>&1
             if [ "$?" -ne 0 ]
              then
                  echo "## ERROR: platypus"
                  exit
             fi
  
             ###samtools
             echo "$(date) Running Sample $SM, Samtools: mpileup & bcftools call"
             /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv samtools mpileup -d8000 -uf ${references["REF"]} -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR ${INPUTFILE} -l ${libraries["REGIONS"]} | bcftools call  -mv -Ov  > $RESULT/$PREFIX\_samtools.vcf  2> $LOG/$PREFIX\_samtools.log
              if [ "$?" -ne 0 ]
                  then
                      echo "## ERROR: samtools mpileup & bcftools call"
                      exit
              fi

             ###freebayes
             echo "$(date) Running Sample $SM, freebayes"
             /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv freebayes  -t ${libraries["REGIONSFB"]} -f ${references["REF"]} -b ${INPUTFILE} | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & DP > 10" -F "LowQual"> $RESULT/$PREFIX\_freebayes.vcf 2> $LOG/$PREFIX\_freebayes.log

             /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  grep '^#' $RESULT/$PREFIX\_freebayes.vcf > $RESULT/$PREFIX\_freebayessorted.vcf && grep -v -e '^#\|^X\|^Y' $RESULT/$PREFIX\_freebayes.vcf | LC_ALL=C sort -t $'\t' -k1,1n -k2,2n >> $RESULT/$PREFIX\_freebayessorted.vcf && grep -e '^X\|^Y' $RESULT/$PREFIX\_freebayes.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> $RESULT/$PREFIX\_freebayessorted.vcf
  
             if [ "$?" -ne 0 ]
                  then
                      echo "## ERROR: freebayes"
                      exit
             fi
             ###pindel

            ###mergevcf
             echo "$(date) Running Sample $SM, Merging"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["GATK"]}/GenomeAnalysisTK.jar -T CombineVariants -R ${references["REF"]} -nt ${params["NT"]} --genotypemergeoption PRIORITIZE -priority haplotypeCaller,platypus,samtools,freebayes  -L ${libraries["REGIONS"]} -V:platypus $RESULT/$PREFIX\_platypus.vcf -V:samtools $RESULT/$PREFIX\_samtools.vcf -V:haplotypeCaller $RESULT/$PREFIX\_haplotypefinal.vcf -V:freebayes $RESULT/$PREFIX\_freebayessorted.vcf -o ${OUTPUTFILE} > $LOG/$PREFIX\_combineVariantCallers.log 2>&1

            INPUTFILE=${OUTPUTFILE}
    }

function functionalAnnotation {
  # lista de archivos vcf separadas por espacios. eg. s1.vcf s2.vcf s3.vcf , va sin comillas ""
  arr=( "${INPUTFILE}")
  for filepath in ${arr[@]}; do
      AUX=$filepath
      filename=${filepath##*/}
      name=${filename%%\.vcf} 
      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 0: Split Alternatives with bptools"
      OUTPUTFILE=$RESULT/$name\_Step0_splitMAA.vcf
    
      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar  -mma  $AUX  $OUTPUTFILE  > $LOG/$name\_Step0_splitMAA.error.log 2>&1

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 1: Annotate with snpEff"
      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step1_SnpEff.vcf
      
      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/snpEff.jar ${references["GRCH_VERSION"]}  -hgvs -t -lof -noStats  -canon -onlyProtein -c ${tools["SNPEFF"]}/snpEff.config $AUX > $OUTPUTFILE 2> $LOG/$name\_Step1_SnpEff.error.log

      if [ "$?" -ne 0 ]
       then
        echo "## ERROR: snpEff command"
        exit
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 2: Annotate with B_Platform VariantAnnotator - select only the highest-impact effect for each variant"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step2_VariantAnnotator.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar -bpann  $AUX  $OUTPUTFILE  > $LOG/$name\_Step2_VariantAnnotator.error.log 2>&1

      if [ "$?" -ne 0 ]
       then
        echo "## ERROR:  Annotate with B_Platform VariantAnnotator"
        exit
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 3: Annotate with dbSNP150"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step3_dbSNP150.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["DBSNP"]} $AUX > $OUTPUTFILE 2> $LOG/$name\_Step3_dbSNP150.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR: Annotate with dbSNP150"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
           OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 4: Annotate with 1000Genomes"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step4_1000Genomes.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["KOMNI"]}  $AUX > $OUTPUTFILE  2> $LOG/$name\_Step4_1000Genomes.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with 1000Genomes"
          read -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi


      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 5: Annotate with hapmap"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step5_hapmap.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["HAPMAP"]}  $AUX > $OUTPUTFILE 2> $LOG/$name\_Step5_hapmap.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with hapmap"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi


      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 6: Annotate with GWASCat (with extra fields)"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step6_GWASCat.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar   ${tools["SNPEFF"]}/SnpSift.jar gwasCat -db ${references["GWASCAT"]} $AUX > $OUTPUTFILE 2> $LOG/$name\_Step6_GWASCat.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with GWASCat (with extra fields)"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
           OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 7: Annotate with dbNSFP"
      #âYou can now specify which fields you want to use for annotation using the '-f' command line option followed by a comma separated list of field names. Defaults fileds are ... â #http://snpeff.sourceforge.net/SnpSift.html#dbNSFP

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step7_dbNSFP.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar dbnsfp -v -f "hg18_pos(1-coor),hg38_chr,hg38_pos,genename,Uniprot_acc,Uniprot_id,cds_strand,CADD_raw,CADD_phred,CADD_raw_rankscore,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_rankscore,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,VEST3_score,VEST3_rankscore,GERP++_RS,GERP++_RS_rankscore,COSMIC_ID,COSMIC_CNT,phyloP46way_placental,phyloP100way_vertebrate,phastCons46way_placental,phastCons100way_vertebrate,Ensembl_geneid,Ensembl_transcriptid"  -db ${references["DBNSFP"]}  $AUX > $OUTPUTFILE 2> $LOG/$name\_Step7_dbNSFP.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with dbNSFP"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 8: Annotate with VarType"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step8_VarType.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o  $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar varType $AUX > $OUTPUTFILE 2>  $LOG/$name\_Step8_VarType.error.log

      if [ "$?" -ne 0 ]
       then
        echo "## ERROR:  Annotate with VarType"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
          then
              exit
          else
              OUTPUTFILE=$AUX
        fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 9: Annotate with Exome Variant Server"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step9_EVS.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate -v  -info MAF ${references["EVS"]} $AUX > $OUTPUTFILE 2> $LOG/$name\_Step9_EVS.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with Exome Variant Server"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
           OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 10: Annotate with PhastCons"
      #http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step10_PhastCons.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar phastCons ${tools["SNPEFF"]}/data/phastCons $AUX >$OUTPUTFILE 2> $LOG/$name\_Step10_PhastCons.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with PhastCons"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 11: Annotate with CADD (Combined Annotation Dependent Depletion) - 1000 Genome variants"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step11_CADD.vcf

     /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate -v ${references["CADD"]} $AUX > $OUTPUTFILE  2> $LOG/$name\_Step11_CADD.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with CADD"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 12: Annotate with ClinVar"

      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step12_clinVar.vcf

      /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info CLNHGVS,CLNALLE,CLNSRC,CLNORIGIN,CLNSRCID,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNACC ${references["CLNVAR"]} $AUX > $OUTPUTFILE  2> $LOG/$name\_Step12_clinVar.error.log

      if [ "$?" -ne 0 ]
       then
          echo "## ERROR:  Annotate with ClinVar"
          read  -p "Do you want to continue? [y,n]" answer
          if [ "$answer" = "n" ]
          then
             exit
          else
             OUTPUTFILE=$AUX
          fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 13: Annotate with PharmGKB"
      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_Step13_pharmGKB.vcf
      
      /usr/bin/time -a -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info PGKB_INDEX,PGKB_GENE,PGKB_DRUG,PGKB_TYPE,PGKB_EVIDENCE,PGKB_DISEASE,PGKB_RACE ${references["PHARMGKB"]} $AUX > $OUTPUTFILE 2> $LOG/$name\_Step13_pharmGKB.error.log

      if [ "$?" -ne 0 ]
       then
         echo "## ERROR:  Annotate with PharmGKB"
         read  -p "Do you want to continue? [y,n]" answer
         if [ "$answer" = "n" ]
         then
            exit
         else
            OUTPUTFILE=$AUX
         fi
      fi

      echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 14: Annotate with ExAC"
      AUX=$OUTPUTFILE
      OUTPUTFILE=$RESULT/$name\_final_annot.vcf
      /usr/bin/time -a -o $TIMELOG/$name\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info AN_Adj,AC_Adj,AC_Het,AC_Hom,AC_Hemi,POPMAX,VQSLOD,GQ_MEAN,GQ_STDDEV,HWP ${references["ExAC"]} $AUX > $OUTPUTFILE 2> $LOG/$name\_Step14_ExAC.error.log

      if [ "$?" -ne 0 ]
       then
         echo "## ERROR:  Annotate with ExAC"
         read  -p "Do you want to continue? [y,n]" answer
         if [ "$answer" = "n" ]
         then
            exit
         else
            OUTPUTFILE=$AUX
         fi
      fi
      INPUTFILE=$OUTPUTFILE
  done
}


function coverage {
    echo -n "Please enter bed file (absolute path from home):"
    read BEDFILE

    echo -n "Please enter bam file (absolute path from home):"
    read BAMFILE
      
    echo -n "Please enter results directory (full path /home/bitgenia/...):"
    read RESULTPATH

    if [[ $BAMFILE == *.bam ]]
     then
      echo -e "#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean" > $RESULTPATH/depthAnalysis.csv
      for region in $(cat "${BEDFILE}" | sed 's/\t/:/' | sed 's/\t/-/'); do
        echo -n "$region" | tr :- \\t
        samtools depth "${BAMFILE}" -a -r $region | awk 'BEGIN {min=100000000; max=0;start=0;s=0};{c++;s+=$3;if(start==0)start=$2;if($3>max)max=$3;if($3<=min)min=$3;l=$2-start}END {print "\011",((l>0)?l:0),"\011",((min<100000000)?min:0),"\011",max,"\011",s/((c>0)?c:1)}'
      done  >> $RESULTPATH/depthAnalysis.csv;
            echo -e "#chrom\tstart\tend\t#readsQueCubrenA\t#basesCubiertasDeA\tsizeOfA\t%OfACubiertas" > $RESULTPATH/breadthAnalysis.csv
            bedtools coverage -a $BEDFILE -b $BAMFILE >> $RESULTPATH/breadthAnalysis.csv;
     else
      echo -e "#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean" >> $RESULTPATH/depthAnalysis.csv
      FILES="$(< $BAMFILE)" 
      for NAME in $FILES; do  
        echo -e "#chrom\tstart\tend\tlength\tmincov\tmaxcov\tmean\t$NAME" >> $RESULTPATH/depthAnalysis.csv  
        for region in $(cat "${BEDFILE}" | sed 's/\t/:/' | sed 's/\t/-/'); do
          echo -n "$region" | tr :- \\t     
          samtools depth "$NAME" -a -r $region | awk 'BEGIN {min=100000000; max=0;start=0;s=0};{c++;s+=$3;if(start==0)start=$2;if($3>max)max=$3;if($3<=min)min=$3;l=$2-start}END {print "\011",((l>0)?l:0),"\011",((min<100000000)?min:0),"\011",max,"\011",s/((c>0)?c:1)}'
        done  >> $RESULTPATH/depthAnalysis.csv; 
      done
        
            FILES="$(< $BAMFILE)"
            for NAME in $FILES; do
               echo -e "$NAME #chrom\tstart\tend\t#readsQueCubrenA\t#basesCubiertasDeA\tsizeOfA\t%OfACubiertas" >> $RESULTPATH/breadthAnalysis.csv    
               bedtools coverage -a $BEDFILE -b "$NAME"
      done >> $RESULTPATH/breadthAnalysis.csv
    fi
}


function generatePathogenicVcf {
  echo -n "Please enter list of genes (space separated):"
  read INPUTFILE

  echo -n "Please enter the file output name:"
  read OUTPUTNAME

  echo -n "Please enter results directory (full path /home/bitgenia/...):"
  read RESULTPATH

  genes=( "${INPUTFILE}")
  pattern="="

 for gene in ${genes[@]}; do
    pattern=$pattern$gene":\|="
  done
  pattern=${pattern::-3}
  echo "^#\|"$pattern
  zgrep "^#\|CLNSIG=[5\|4]\||5\||4" ${GRCh37["CLNVAR"]} | grep "^#\|"$pattern > $RESULTPATH/$OUTPUTNAME.vcf

}

function vcf2bed {
  echo -n "Please enter vcf file (absolute path from home):"
  read INPUTFILE
  
  INPUTFILE=$(trim ${INPUTFILE})
  if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
  then
    echo "Error: vcf file is null or not exists"
    exit
  fi
  
  echo -n "Please enter results directory (full path /home/bitgenia/...):"
  read RESULTPATH

  RESULTPATH=$(trim ${RESULTPATH})
  if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
   then
      echo "Error: the result directory is null or not exists"
      exit
  fi

  grep -v "^#" "${INPUTFILE}" | awk '{if ($2 > 1000) print $1, $2-1000, $2+1000;  else print $1, 0, $2+1000; }' > $RESULTPATH/vcfVariants.bed
}

function reduceBam {
  
  echo -n "Please enter the sample name:"
  read SAMPLENAME

  SAMPLENAME=$(trim ${SAMPLENAME})
  if [ -z "$SAMPLENAME" ]
  then
    echo "Error: Sample name is null or empty"
    exit
  fi

  echo -n "Please enter bam file (absolute path from home):"
  read INPUTFILE

  INPUTFILE=$(trim ${INPUTFILE})
  if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
  then
    echo "Error: bam file is null or not exists"
    exit
  fi
  
  echo -n "Please enter bed file (absolute path from home):"
  read INPUTFILEBED

  INPUTFILEBED=$(trim ${INPUTFILEBED})
  if [ -z "$INPUTFILEBED" ] || [ ! -f "$INPUTFILEBED" ]
  then
    echo "Error: bed file is null or not exists"
    exit
  fi

  echo -n "Please enter results directory (full path /home/bitgenia/...):"
  read RESULTPATH

  RESULTPATH=$(trim ${RESULTPATH})
  if [ -z "$RESULTPATH" ] || [ ! -d "$RESULTPATH" ]
   then
      echo "Error: the result directory is null or not exists"
      exit
  fi

  echo samtools view -bh -L "${INPUTFILEBED}" "${INPUTFILE}" > "${RESULTPATH}"/"${SAMPLENAME}"_reduced.bam
  samtools view -bh -L "${INPUTFILEBED}" "${INPUTFILE}" > "${RESULTPATH}"/"${SAMPLENAME}"_reduced.bam
  samtools index "${RESULTPATH}"/"${SAMPLENAME}"_reduced.bam
}

function variantDiscoveryBPResolution {

    echo "$(date) Running Sample $SM, GATK: HaplotypeCaller in gVCF mode"
    INPUTFILEK=${INPUTFILE}
    OUTPUTFILEK=$RESULT/$PREFIX\_haplotypecaller.g.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T HaplotypeCaller -L ${references["KNOWNVARS"]} -nct ${params["NCT"]} -I ${INPUTFILEK} -o ${OUTPUTFILEK} --emitRefConfidence BP_RESOLUTION  > $LOG/$PREFIX.haplotypecallerGVCF.log 2>&1

    echo "$(date) Running Sample $SM, GATK: Genotype gVCF"
    #Get a vcf from a gvcf
    INPUTFILEK=${OUTPUTFILEK}
    OUTPUTFILEK=$RESULT/$PREFIX\_genotypegvcf.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T GenotypeGVCFs -G StandardAnnotation -o ${OUTPUTFILEK} -allSites --variant ${INPUTFILEK} > $LOG/$PREFIX.genotypeGVCF.log 2>&1
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T GenotypeGVCFs -G StandardAnnotation -o ${OUTPUTFILEK} -L ${references["KNOWNVARS"]}  --variant ${INPUTFILEK} > $LOG/$PREFIX.genotypeGVCF.log 2>&1


    echo "$(date) Running Sample $SM, Split annotate with bptools: FORMAT#RGQ<30"
    INPUTFILEK=${OUTPUTFILEK}
    OUTPUTFILEK=$RESULT/$PREFIX\_rgq_known_variant_calling.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar  -annotate -exp "FORMAT#RGQ<30" -tag "FILTER#LowConfidence" ${INPUTFILEK} ${OUTPUTFILEK}  > $LOG/$PREFIX\_annotate_error.log 2>&1


    echo "$(date) Running Sample $SM, Split annotate with bptools: FORMAT#DP<10"
    INPUTFILEK=${OUTPUTFILEK}
    OUTPUTFILEK=$RESULT/$PREFIX\_dp_known_variant_calling.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar  -annotate -exp "FORMAT#DP<10" -tag "FILTER#LowConfidence" ${INPUTFILEK} ${OUTPUTFILEK}  > $LOG/$PREFIX\_annotate_error.log 2>&1

    echo "$(date) Running Sample $SM, Split annotate with bptools: FORMAT#GT=./."
    INPUTFILEK=${OUTPUTFILEK}
    OUTPUTFILEK=$RESULT/$PREFIX\_final_known_variant_calling.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar  -annotate -exp "FORMAT#GT=./." -tag "FILTER#MissingData" ${INPUTFILEK}  ${OUTPUTFILEK} > $LOG/$PREFIX\_annotate_error.log 2>&1

}

function sendToFTP {

 if [ ! -f "$INPUTFILE" ]
     then
        #lista de archivos vcf separadas por espacios. eg. s1.vcf s2.vcf s3.vcf , va sin comillas ""
        echo -n "Please enter vcf file to upload (absolute path from home):"
        read INPUTFILE
        arr=( "${INPUTFILE}")
        for filepath in ${arr[@]}; do
          echo "$filepath"
          filepath=$(trim ${filepath})
          if [ -z "$filepath" ] || [ ! -f "$filepath" ]
          then
            echo "Error: $filepath file is null or not exists"
            exit
          fi
        done
 fi

 if [ -z "$DOMAIN" ]
     then
       echo  -n "Enter code domain:"
       read DOMAIN
       DOMAIN=$(trim ${DOMAIN})
 fi

 if [ -n "${DOMAIN}" ]&&[ "${DOMAIN}" -ne 0 ]
    then
      echo "####  uploading vcf to Bplatform ftp #######"
      #zip -j $INPUTFILE\.zip  $INPUTFILE
      sshpass -p "${PASS}" scp  $INPUTFILE "${USR}"@"${SERVER}":/home/"${USR}"/upload/"${DOMAIN}";
      if [ "$?" -ne 0 ]
       then
           echo "Upload Unsuccessful!"
           exit
      fi
  fi
}


function resetVariables {
    INPUTFILE="";
    OUTPUTFILE="";
    DOMAIN="";
}

function interactiveScript {

    while(true)
    do
        echo
        printf "Choose from the following operations: \n"
        printf "[p]reprocessing\n"
        printf "[v]ariantDiscovery\n"
        printf "[m]erge\n"
        printf "[f]unctionalAnnotation\n"
        printf "[l]ClinVariants\n"
        printf "[k]nownvariant\n"
        printf "[c]ompletePipeline\n"
        printf "[u]pload file to Bplatform\n"
        read -p "Your choice: " choice
        echo   "-----------------------------------------"
        case $choice in
        [pP])
            getReferences
            getInteractiveSampleInfo
            preprocessing
        ;;
        [vV])
            getReferences
            getInteractiveBamInfo
            variantDiscovery
        ;;
        [mM])
            getReferences
            getInteractiveBamInfo
            fullvariantDiscovery
        ;;
        [fF])
            getReferences
            getInteractiveVcfInfo
            functionalAnnotation
        ;;
        [lL])
           generatePathogenicVcf
        ;;
        [kK])
            getReferences
            getInteractiveBamInfo "k"
            variantDiscoveryBPResolution
        ;;
        [cC])
            getReferences
            getInteractiveSampleInfo
            preprocessing
            variantDiscovery "h"
            functionalAnnotation
            sendToFTP
        ;;
        [uU])
            sendToFTP
        ;;
        [bB])
            vcf2bed
        ;;
        [rR])
            reduceBam
        ;;
        [dD])
            coverage
        ;;
        *)
           echo "wrong choice!"
           exit;
        esac
        resetVariables
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

function displayHelp {

    printf "usage: exome [OPTION] [SAMPLENAME] [fastq1 .. fastq2] [outdir] [domain]\n"
    printf "parameters:\n"
    printf "            h --help       Prints this help dialog.\n"
    printf "            d --default    Non-intereactive shell script\n"
    echo   "--------------------------------------------------------"
    printf "examples:\n"
    echo
    printf "Non-interactive mode\n"
    echo
    printf "    bash shell.bash d 888 /home/user/exome/888_1.fastq /home/user/exome/888_2.fastq /home/user/exome/888_results\n"
    printf "    bash shell.bash d 999 /home/user/exome/999_1.fastq /home/user/exome/999_2.fastq /home/user/exome/999_results  99 \n"
    echo
    printf "Interactive mode\n"
    echo
    printf "    bash shell.bash \n"
    echo

}
########################### MAIN ######################
 echo
 echo   "---------------------------------------------"
 echo   "   Welcome to Analysis Pipeline ${VERSION} - ${UPDATE}"
 echo   "   For support and documentation contact: info@bitgenia.com"
 echo   "   www.bitgenia.com - www.biargentina.com.ar"
 echo   "   Copyright (c) 2013-2016 BITGENIA - BIA"
 echo   "---------------------------------------------"

OPTION=$(trim ${1})
if [ -z "${OPTION}" ]
   then
       interactiveScript
   else
       case "${OPTION}" in
        [dD])
          getNonInteractiveSampleInfo ${2} ${3} ${4} ${5} ${6} ${7}
          processLanes
#          preprocessing
          variantDiscovery "h"
          functionalAnnotation
          sendToFTP
        ;;
        [mM])
          getNonInteractiveSampleInfo ${2} ${3} ${4} ${5} ${6}
          preprocessing
          fullvariantDiscovery
          functionalAnnotation
          sendToFTP
        ;;
        [hH])
          displayHelp
        ;;
        *)
           echo "wrong choice!"
           exit
        esac
fi
