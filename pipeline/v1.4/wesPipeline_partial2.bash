#!/bin/bash

## Script version
VERSION="1.4"

## DECLARE VARIABLES
USR='sample'
PASS='sample321654'
SERVER=192.168.210.77

declare -A tools=(["BWA"]="/home/bitgenia/tools/aligners/bwa-0.7.12" ["BT2"]="/home/bitgenia/tools/aligners/bowtie2-2.2.6" ["BBMAP"]="/home/bitgenia/tools/aligners/bbmap" ["PICARD"]="/home/bitgenia/tools/accessoryTools/picard-tools-2.0.1" ["GATK"]="/home/bitgenia/tools/variantCallers" ["PRINSEQ"]="/home/bitgenia/tools/accessoryTools/prinseq-lite-0.20.4" ["PLATYPUS"]="/home/bitgenia/tools/variantCallers/Platypus_0.8.1" ["SNPEFF"]="/home/bitgenia/tools/functionalAnn/snpEff" ["BPLATANN"]="/home/bitgenia/tools/accessoryTools" ["FASTQC"]="/home/bitgenia/tools/accessoryTools/FastQC" ["SAMTOOLS"]="/home/bitgenia/tools/accessoryTools/samtools-1.3");


declare -A GRCh37=(["HGREF"]="/home/bitgenia/bundle/hgref" ["BT2_REF"]="/home/bitgenia/bundle/hgref/human_g1k_v37_decoy" ["REF"]="/home/bitgenia/bundle/hgref/human_g1k_v37_decoy.fasta" ["KGINDELS"]="/home/bitgenia/bundle/hgref/1000G_phase1.indels.b37.vcf" ["MILLSINDELS"]="/home/bitgenia/bundle/hgref/Mills_and_1000G_gold_standard.indels.b37.vcf" ["DBSNP"]="/home/bitgenia/bundle/hgref/dbSNP_144_GRCH37.vcf.gz" ["GWASCAT"]="/home/bitgenia/bundle/hgref/gwascatalog.txt"  ["GWASCAT_ADJUSTED"]="/home/bitgenia/bundle/hgref/gwascatalog_adjusted.vcf" ["EVS"]="/home/bitgenia/bundle/hgref/ESP6500SI-V2-SSA137.snps_indels.vcf" ["CADD"]="/home/bitgenia/bundle/hgref/CADD_1000G.vcf" ["CLNVAR"]="/home/bitgenia/bundle/hgref/clinvar_20160502.vcf" ["PHARMGKB"]="/home/bitgenia/bundle/hgref/PharmGKBvcf_2016.vcf" ["ExAC"]="/home/bitgenia/bundle/hgref/ExAC.r0.3.1.sites.vep.vcf.gz" ["HAPMAP"]="/home/bitgenia/bundle/hgref/hapmap_3.3.b37.vcf" ["KOMNI"]="/home/bitgenia/bundle/hgref/1000G_omni2.5.b37.vcf" ["KNOWNVARS"]="/home/bitgenia/scripts/known_variants_analysis-devel/grch37mybaby525.vcf" ["DBNSFP"]="/home/bitgenia/bundle/hgref/dbNSFP2.9.txt.gz" ["GRCH_VERSION"]="GRCh37.75");


declare -A hg19=(["HGREF"]="/home/bitgenia/bundle/hg19" ["BT2_REF"]="/home/bitgenia/bundle/hg19" ["REF"]="/home/bitgenia/bundle/hg19/ucsc.hg19.fasta" ["KGINDELS"]="/home/bitgenia/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf" ["MILLSINDELS"]="/home/bitgenia/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf" ["DBSNP"]="/home/bitgenia/bundle/hg19/dbsnp_144.hg19.vcf.gz" ["GWASCAT"]="/home/bitgenia/bundle/hg19/gwascatalog.txt" ["GWASCAT_ADJUSTED"]="/home/bitgenia/bundle/hg19/gwascatalog_adjusted.vcf" ["EVS"]="/home/bitgenia/bundle/hg19/ESP6500SI-V2-SSA137.snps_indels.vcf" ["CADD"]="/home/bitgenia/bundle/hg19/CADD_1000G.vcf" ["CLNVAR"]="/home/bitgenia/bundle/hg19/clinvar_20160502.vcf" ["PHARMGKB"]="/home/bitgenia/bundle/hg19/PharmGKBvcfhg19_2016.vcf" ["ExAC"]="/home/bitgenia/bundle/hg19/ExAC.r0.3.1.sites.vep.vcf.gz" ["HAPMAP"]="/home/bitgenia/bundle/hg19/hapmap_3.3.hg19.sites.vcf" ["KOMNI"]="/home/bitgenia/bundle/hg19/1000G_omni2.5.hg19.sites.vcf" ["KNOWNVARS"]="/home/bitgenia/scripts/known_variants_analysis-devel/clinvar_Result/hg19mybaby525.vcf" ["DBNSFP"]="/home/bitgenia/bundle/hg19/dbNSFP2.9.txt.gz" ["GRCH_VERSION"]="hg19");
#XMX=Java heap max size
#NT= Num of data threads sent to processor, Http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
#NCT= Num of CPU threads for each data thread, http://gatkforums.broadinstitute.org/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
#QFILTER_READS PARAMETER OPTIONAL min.qual.mean: "Minimum mean quality" TYPE INTEGER(Filter reads with quality score mean below the given value.)
#MINPRUN minPruning sets the kmer size for the BRUIJN tree in variant calling.
#CALL Variant Call Confidence when creating .vcf file.
#EMIT  Variant Emition Confidence when creating .vcf file.
declare -A params=(["XMX"]="16g" ["NT"]="11" ["NCT"]="8" ["QFILTER_READS"]="10" ["MINPRUN"]="3" ["CALL"]="30.0" ["EMIT"]="10.0");


#LB= RGLB=String   REQUIRED Library.
#PL=     # RG PL=String  REQUIRED platform (e.g. illumina, solid). IMPORTANTE para GATK: no puede ponerse UNKNOWN por más que lo muestre como posible. Debe especificarse una plataforma, la marca en general (ILLUMINA y no HISEQ200, por ejemplo).
#PU=     # RGPU=String   REQUIRED platform unit (eg. run barcode).
#ID=1    # RGID=String   ID Default value: 1. This option can be set to 'null' to clear the default value.
#CN=null # RGCN=String   sequencing center name Default value: null.
#DS=null # RGDS=String   description Default value: null.
#DT=null # RGDT=Iso8601Date      run date Default value: null.
#PI=null # RGPI=INTEGER  predicted insert size Default value: null.
declare -A readGroups=(["LB"]="1" ["PL"]="ILLUMINA" ["PU"]="unknown" ["ID"]="1" ["CN"]="" ["DS"]="" ["DT"]="" ["PI"]="");




#[design ID]_Padded.bed - This BED file contains a single track of the genomic regions that you can expect to sequence when using the design for target enrichment. To determine these regions, the program extends the regions in the Covered BED file by 100 bp on each side.
declare -A librarySSV4UTRsGrch37=(["REGIONS"]="/home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV4+UTRs_Padded.bed" ["REGIONSTXT"]="");
declare -A librarySSV5UTRsGrch37=(["REGIONS"]="/home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV5+UTRs_Padded.bed" ["REGIONSTXT"]="");
declare -A librarySSV5Grch37=(["REGIONS"]="/home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV5_Padded.bed" ["REGIONSTXT"]="/home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV5_Regions.txt" ["REGIONSFB"]="/home/bitgenia/bundle/hgref/library/SureSelectHumanAllExonV5_Targets.bed");

declare -A libraryTruSightOneGrch37=(["REGIONS"]="/home/bitgenia/bundle/hgref/library/TruSight_One_v1.1.bed" ["REGIONSTXT"]="");


declare -A librarySSV4UTRshg19=(["REGIONS"]="/home/bitgenia/bundle/hg19/library/SureSelectHumanAllExonV4+UTRs_Padded.bed" ["REGIONSTXT"]="" ["REGIONSFB"]="");

declare -A librarySSV5UTRshg19=(["REGIONS"]="/home/bitgenia/bundle/hg19/library/SureSelectHumanAllExonV5+UTRs_Padded.bed" ["REGIONSTXT"]="" ["REGIONSFB"]="");

declare -A librarySSV5hg19=(["REGIONS"]="/home/bitgenia/bundle/hg19/library/SureSelectHumanAllExonV5_Padded.bed" ["REGIONSTXT"]="/home/bitgenia/bundle/hg19/library/SureSelectHumanAllExonV5_Regions.txt" ["REGIONSFB"]="/home/bitgenia/bundle/hg19/library/SureSelectHumanAllExonV5_Targets.bed");

declare -A libraryTruSightOnehg19=(["REGIONS"]="/home/bitgenia/bundle/hg19/library/TruSight_One_v1.1.bed" ["REGIONSTXT"]="");

declare -A libraries=();

declare -A references=();

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
            printf "[a]GRch37\n"
            printf "[b]hg19\n"
            read -p "Your choice: " choice
            echo   "-----------------------------------------"
            echo
            printf "Choose from the following capture prep Kit:\n"
            printf "[a]SSV5+UTRs\n"
            printf "[b]SSV4+UTRs\n"
            printf "[c]SSV5\n"
            printf "[e]TruSight One"

            read -p "Your choice: " choiceK
            echo
      fi
      case $choice in
          [aA])
            for K in "${!GRCh37[@]}"; do references[$K]=${GRCh37[$K]}; done
            case $choiceK in
              [a])
                for K in "${!librarySSV5UTRsGrch37[@]}"; do libraries[$K]=${librarySSV5UTRsGrch37[$K]}; done
              ;;
              [b])
                for K in "${!librarySSV4UTRsGrch37[@]}"; do libraries[$K]=${librarySSV4UTRsGrch37[$K]}; done
              ;;
              [c])
               for K in "${!librarySSV5Grch37[@]}"; do libraries[$K]=${librarySSV5Grch37[$K]}; done
              ;;
              [e])
               for K in "${!libraryTruSightOneGrch37[@]}"; do libraries[$K]=${libraryTruSightOneGrch37[$K]}; done
              ;;
              *)
                 echo "wrong choice!"
                 exit;
            esac
          ;;
          [bB])
            for K in "${!hg19[@]}"; do references[$K]=${hg19[$K]}; done
            case $choiceK in
              [a])
                for K in "${!librarySSV5UTRshg19[@]}"; do libraries[$K]=${librarySSV5UTRshg19[$K]}; done
              ;;
              [b])
                for K in "${!librarySSV4UTRshg19[@]}"; do libraries[$K]=${librarySSV4UTRshg19[$K]}; done
              ;;
              [c])
                for K in "${!librarySSV5hg19[@]}"; do libraries[$K]=${librarySSV5hg19[$K]}; done
              ;;
              [e])
               for K in "${!libraryTruSightOnehg19[@]}"; do libraries[$K]=${libraryTruSightOnehg19[$K]}; done
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


function getNonInteractiveSampleInfo {

    getReferences "a" "e"
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

    DOMAIN=$(trim ${5})
    if [[ -z "${DOMAIN}" ]]
     then
        echo "Warning: without domain id we can't not upload vcf file to B_Platform's ftp"
        DOMAIN=0
    fi
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


function preprocessing {

    #echo "$(date) Running Sample $SM, FastQC"
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv ${tools["FASTQC"]}/fastqc -t 6 $FIRSTREADS $MATES --outdir=$RESULT  > $LOG/$PREFIX\_fastqc.log 2>&1

    #if [ "$?" -ne 0 ]
    # then
    #  echo "## ERROR: FastQC command"
    #  exit
    #fi

    #echo "$(date) Running Sample $SM, PrinSEQ-lite"
    #Quality control. Input: 2 fastq files (one per mate); output: 2 fastq files corrected with qc methods
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  perl ${tools["PRINSEQ"]}/prinseq-lite.pl -fastq $FIRSTREADS -fastq2 $MATES -out_good $RESULT/$PREFIX_QC -out_bad $RESULT/$PREFIX_QC\_BAD -min_qual_mean ${params["QFILTER_READS"]} -log $LOG/$PREFIX\_prinseq.log > $LOG/$PREFIX\_prinseq_results.log 2>&1

    #if [ "$?" -ne 0 ]
    # then
    #  echo "## ERROR: PrinSEQ command"
    #  exit
    #fi

    #echo "$(date) Running Sample $SM, BWA: mem"
    #OUTPUTFILE=$RESULT/$PREFIX\_aligned.sam
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv ${tools["BWA"]}/bwa mem -M -R "@RG\tID:${readGroups[ID]}\tSM:$SM\tLB:${readGroups[LB]}\tPL:${readGroups[PL]}" -t ${params["NT"]} ${references["REF"]} $RESULT/$PREFIX_QC\_1.fastq $RESULT/$PREFIX_QC\_2.fastq > ${OUTPUTFILE} 2> $LOG/$PREFIX\_bwa_mem.log

    #if [ "$?" -ne 0 ]
    # then
    #  echo "## ERROR: BWA command"
    #  exit
    #fi

    #echo "$(date) Running Sample $SM, PICARD: CleanSam"
    #INPUTFILE=${OUTPUTFILE}
    #OUTPUTFILE=$RESULT/$PREFIX\_cleaned.sam
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar CleanSam I=${INPUTFILE} O=${OUTPUTFILE} CREATE_INDEX=true > $LOG/$PREFIX.cleansam.log 2>&1
    #if [ "$?" -ne 0 ]
    # then
    #   echo "## ERROR: Picard CleanSam"
    #  exit
    #fi

    #echo "$(date) Running Sample $SM, PICARD: SortSam"
    #INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX\_sorted.bam
    #/usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar SortSam I=${INPUTFILE} O=${OUTPUTFILE} SO=coordinate CREATE_INDEX=true > $LOG/$PREFIX.sortsam.log 2>&1

    #if [ "$?" -ne 0 ]
    # then
    #   echo "## ERROR: Picard SortSam"
    #  exit
    #fi

    echo "$(date) Running Sample $SM, SamTools: FlagStat"
    ${tools["SAMTOOLS"]}/samtools flagstat ${OUTPUTFILE} > $RESULT/$PREFIX\_PreMarkDuplicatesStats.txt 2>&1

    echo "$(date) Running Sample $SM, PICARD: MarkDuplicates"
    INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX\_marked.bam
    #Marcar duplicados. Input: bam, output: bam.
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true AS=true MAX_RECORDS_IN_RAM=10000000 I=${INPUTFILE} O=${OUTPUTFILE} M=$RESULT/$PREFIX\_dedupped.metrics.txt CREATE_INDEX=true > $LOG/$PREFIX.markduplicates1.log 2>&1

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard markduplicates"
      exit
    fi

    echo "$(date) Running Sample $SM, SamTools: FlagStat"
    ${tools["SAMTOOLS"]}/samtools flagstat ${OUTPUTFILE} > $RESULT/$PREFIX\_PostMarkDuplicatesStats.txt 2>&1

    echo "$(date) Running Sample $SM, PICARD: CollectAlignmentSummaryMetrics"
    INPUTFILE=${OUTPUTFILE}
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["PICARD"]}/picard.jar CollectAlignmentSummaryMetrics I=${INPUTFILE} O=$RESULT/$PREFIX\_collectAlignmentSummaryMetrics.txt R=${references["REF"]} CREATE_INDEX=true > $LOG/$PREFIX.collectAlignmentSummaryMetrics.log 2>&1

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR: Picard CollectAlignmentSummaryMetrics"
      exit
    fi

    #In order to reduce the number of miscalls of INDELs in your data it is helpful to realign your raw gapped alignment with the Broad’s GATK (https://www.broadinstitute.org/gatk/) Realigner.

    echo "$(date) Running Sample $SM, GATK: RealignerTargetCreator"
    INPUTFILE=${OUTPUTFILE}
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${references["REF"]} -nt ${params["NT"]} --known ${references["KGINDELS"]} --known ${references["MILLSINDELS"]} -I ${INPUTFILE} -o $RESULT/$PREFIX.intervals > $LOG/$PREFIX.realignertargetcreator.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: GATK: RealignerTargetCreator"
      exit
    fi

    echo "$(date) Running Sample $SM, GATK: IndelRealigner"
    INPUTFILE=${OUTPUTFILE}
    OUTPUTFILE=$RESULT/$PREFIX\_realigned.bam
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T IndelRealigner -R ${references["REF"]} -known ${references["KGINDELS"]} -known ${references["MILLSINDELS"]} -I ${INPUTFILE} -targetIntervals $RESULT/$PREFIX.intervals -o ${OUTPUTFILE} > $LOG/$PREFIX.indelrealigner.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: GATK: IndelRealigner"
      exit
    fi

    #BQSR from the Broad’s GATK allows you to reduce the effects of analysis artefacts produced by your sequencing machines. It does this in two steps, the first analyses your data to detect covariates and the second compensates for those covariates by adjusting quality scores.

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
    OUTPUTFILE=$RESULT/$PREFIX\_bqsr.bam
    #https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_readutils_PrintReads.php
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -T PrintReads -R ${references["REF"]} -I ${INPUTFILE} -BQSR $RESULT/$PREFIX.recal.table -nct ${params["NCT"]}  -o ${OUTPUTFILE} > $LOG/$PREFIX.printreads.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: GATK: PrintReads"
      exit
    fi
  INPUTFILE=$RESULT/$PREFIX\_bqsr.bam
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
        echo -n "Please enter indexed bam file (absolute path from home):"
        read INPUTFILE

        INPUTFILE=$(trim ${INPUTFILE})
        if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
        then
          echo "Error: vcf file is null or not exists"
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
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv ${tools["SAMTOOLS"]}/samtools mpileup -d8000 -uf ${references["REF"]} ${INPUTFILE} -l ${libraries["REGIONS"]} | bcftools call  -mv -Ov  > ${OUTPUTFILE}
            if [ "$?" -ne 0 ]
                then
                    echo "## ERROR: samtools mpileup & bcftools call"
                    exit
            fi
         ;;
        [fF])
            echo "$(date) Running Sample $SM, freebayes"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv freebayes  -t ${libraries["REGIONSFB"]} -f ${references["REF"]} -b ${INPUTFILE} | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & DP > 20" -F "LowQual"> ${OUTPUTFILE}
            if [ "$?" -ne 0 ]
                then
                    echo "## ERROR: freebayes"
                    exit
            fi
         ;;
        [hH])
            echo "$(date) Running Sample $SM, GATK: HaplotypeCaller"
            /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["GATK"]}/GenomeAnalysisTK.jar -R ${references["REF"]} -T HaplotypeCaller -L ${libraries["REGIONS"]} -nct ${params["NCT"]} -I ${INPUTFILE} --dbsnp ${references["DBSNP"]} -o $RESULT/$PREFIX\_haplotypecaller.vcf -stand_call_conf ${params["CALL"]} -stand_emit_conf ${params["EMIT"]} --minPruning ${params["MINPRUN"]} > $LOG/$PREFIX.haplotypecaller.log 2>&1
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

            sed "s/\${SM}.variant/\${SM}/" $RESULT/$PREFIX\_removesamplevariant2.vcf > ${OUTPUTFILE} 2> $LOG/$PREFIX\_final2_sed.log
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
        echo -n "Please enter vcf file (absolute path from home):"
        read INPUTFILE

        INPUTFILE=$(trim ${INPUTFILE})
        if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
        then
          echo "Error: vcf file is null or not exists"
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
    fi
}

function functionalAnnotation {

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 0: Split Alternatives with bptools"
    OUTPUTFILE=$RESULT/$PREFIX\_vcf_maa_annot.vcf
    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar  -mma  $INPUTFILE  $OUTPUTFILE  > $LOG/$PREFIX\_Step0_splitMAA.error.log 2>&1

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 1: Annotate with snpEff"
    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_vcf_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/snpEff.jar  ${references["GRCH_VERSION"]} -hgvs -t -lof -noStats -c ${tools["SNPEFF"]}/snpEff.config $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step1_SnpEff.error.log

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR: snpEff command"
      exit
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 2: Annotate with B_Platform VariantAnnotator - select only the highest-impact effect for each variant"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_snpEff_VAnnot_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["BPLATANN"]}/bptools.jar -bpann  $INPUTFILE  $OUTPUTFILE  > $LOG/$PREFIX\_Step2_VariantAnnotatior.error.log 2>&1

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR:  Annotate with B_Platform VariantAnnotator"
      exit
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 3: Annotate with dbSNP147"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["DBSNP"]} $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step3_dbSNP147.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR: Annotate with dbSNP147"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
         OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 4: Annotate with 1000Genomes"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["KOMNI"]}  $INPUTFILE > $OUTPUTFILE  2> $LOG/$PREFIX\_Step4_1000Genomes.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with 1000Genomes"
        read -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi


    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 5: Annotate with hapmap"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate ${references["HAPMAP"]}  $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step5_hapmap.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with hapmap"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi


    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 6: Annotate with GWASCat (with extra fields)"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar   ${tools["SNPEFF"]}/SnpSift.jar gwasCat -v -db ${references["GWASCAT"]} $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step6_GWASCat.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with GWASCat (with extra fields)"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
         OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 7: Annotate with dbNSFP"
    #“You can now specify which fields you want to use for annotation using the '-f' command line option followed by a comma separated list of field names. Defaults fileds are ... “ #http://snpeff.sourceforge.net/SnpSift.html#dbNSFP

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar dbnsfp -v -f "rs_dbSNP141,genename,Uniprot_acc,Uniprot_id,Uniprot_aapos,cds_strand,refcodon,codonpos,Ensembl_geneid,Ensembl_transcriptid,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,clinvar_rs,clinvar_clnsig,clinvar_trait" -db ${references["DBNSFP"]}  $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step7_dbNSFP.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with dbNSFP"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 8: Annotate with VarType"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o  $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar varType $INPUTFILE > $OUTPUTFILE 2>  $LOG/$PREFIX\_Step8_VarType.error.log

    if [ "$?" -ne 0 ]
     then
      echo "## ERROR:  Annotate with VarType"
      read  -p "Do you want to continue? [y,n]" answer
      if [ "$answer" = "n" ]
        then
            exit
        else
            OUTPUTFILE=$INPUTFILE
      fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 9: Annotate with Exome Variant Server"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_EVS_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv  java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate -v  -info MAF ${references["EVS"]} $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step9_EVS.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with Exome Variant Server"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
         OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 10: Annotate with PhastCons"
    #http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_PhastCons_EVS_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar phastCons ${tools["SNPEFF"]}/data/phastCons $INPUTFILE >$OUTPUTFILE 2> $LOG/$PREFIX\_Step10_PhastCons.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with PhastCons"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi



    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 11: Annotate with CADD (Combined Annotation Dependent Depletion) - 1000 Genome variants"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_CADD_PhastCons_EVS_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar  ${tools["SNPEFF"]}/SnpSift.jar annotate -v ${references["CADD"]} $INPUTFILE > $OUTPUTFILE  2> $LOG/$PREFIX\_Step11_CADD.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with CADD"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 12: Annotate with ClinVar"

    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_CLNVAR_CADD_PhastCons_EVS_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf

    /usr/bin/time -a -f "%C,%e,%K,%M,%t,%I/%O,%P,%W,%c,%w" -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info CLNHGVS,CLNALLE,CLNSRC,CLNORIGIN,CLNSRCID,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNACC ${references["CLNVAR"]} $INPUTFILE > $OUTPUTFILE  2> $LOG/$PREFIX\_Step12_clinVar.error.log

    if [ "$?" -ne 0 ]
     then
        echo "## ERROR:  Annotate with ClinVar"
        read  -p "Do you want to continue? [y,n]" answer
        if [ "$answer" = "n" ]
        then
           exit
        else
           OUTPUTFILE=$INPUTFILE
        fi
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 13: Annotate with PharmGKB"
    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_PharmGKB_CLNVAR_CADD_PhastCons_EVS_VarType_dbNSFP_GWASCat_hapmap_1000G_dbSNP147_snpEff_annot.vcf
    /usr/bin/time -a -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info PGKB_INDEX,PGKB_GENE,PGKB_DRUG,PGKB_TYPE,PGKB_EVIDENCE,PGKB_DISEASE,PGKB_RACE ${references["PHARMGKB"]} $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step13_pharmGKB.error.log

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR:  Annotate with PharmGKB"
       exit
      else
       INPUTFILE=$OUTPUTFILE
    fi

    echo "$(date) Running Sample $SM, SnpEff Annotation Script - Step 14: Annotate with ExAC"
    INPUTFILE=$OUTPUTFILE
    OUTPUTFILE=$RESULT/$PREFIX\_final_annot.vcf
    /usr/bin/time -a -o $TIMELOG/$PREFIX\_time.log.csv java -Xmx${params["XMX"]} -jar ${tools["SNPEFF"]}/SnpSift.jar annotate -v -info AN_Adj,AC_Adj,AC_Het,AC_Hom,AC_Hemi,POPMAX,VQSLOD,GQ_MEAN,GQ_STDDEV,HWP ${references["ExAC"]} $INPUTFILE > $OUTPUTFILE 2> $LOG/$PREFIX\_Step14_ExAC.error.log

    if [ "$?" -ne 0 ]
     then
       echo "## ERROR:  Annotate with ExAC"
       exit
      else
       INPUTFILE=$OUTPUTFILE
    fi
   #SnpSift.jar annotate ../bundle/hgref/ExAC.r0.3.1.sites.vep.vcf.gz test.vcf > test-exac.vcf
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
        echo -n "Please enter vcf file to upload (absolute path from home):"
        read INPUTFILE

        INPUTFILE=$(trim ${INPUTFILE})
        if [ -z "$INPUTFILE" ] || [ ! -f "$INPUTFILE" ]
            then
                echo "Error: vcf file is null or not exists"
                exit
        fi
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
      sshpass -p "${PASS}" scp  $INPUTFILE  "${USR}"@"${SERVER}":/home/"${USR}"/upload/"${DOMAIN}";
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
        printf "[f]unctionalAnnotation\n"
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
        [fF])
            getReferences
            getInteractiveVcfInfo
            functionalAnnotation
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
            #variantDiscoveryBPResolution&
            variantDiscovery "h"
            functionalAnnotation
            sendToFTP
        ;;
        [uU])
            sendToFTP
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
    printf "    bash shell.sh d 888 /home/user/exome/888_1.fastq /home/user/exome/888_2.fastq /home/user/exome/888_results\n"
    printf "    bash shell.sh d 999 /home/user/exome/999_1.fastq /home/user/exome/999_2.fastq /home/user/exome/999_results  99 \n"
    echo
    printf "Interactive mode\n"
    echo
    printf "    bash shell.sh \n"
    echo

}
########################### MAIN ######################
 echo
 echo   "---------------------------------------------"
 echo   "   Welcome to Analysis Pipeline ${VERSION}"
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
          getNonInteractiveSampleInfo ${2} ${3} ${4} ${5} ${6}
          preprocessing
          #variantDiscoveryBPResolution&
          variantDiscovery "h"
          functionalAnnotation
          sendToFTP
        ;;
        [hH])
          displayHelp
        ;;
        *)
           echo "wrong choice!"
           exit;
        esac
fi
