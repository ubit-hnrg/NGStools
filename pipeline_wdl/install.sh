#!/bin/bash

pipelineVersion="HNRG-pipeline-V0.1"
#tools
cromwell='https://github.com/broadinstitute/cromwell/releases/download/37/cromwell-37.jar'
wdltool='https://github.com/broadinstitute/wdltool/releases/download/0.14/wdltool-0.14.jar'
bwakit='http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download'
fastp='http://opengene.org/fastp/fastp'

# GATK version
gatkVersion='4.0.6.0'
gatklink='https://github.com/broadinstitute/gatk/releases/download/'$gatkVersion'/gatk-'$gatkVersion'.zip'

#minimal required dbs
#indels Mills&1000G
Mils_1000G_b37_vcf='ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz'                # required for preprocessing
#dbsnp b37
dbsnp_b37_vcf='ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz'



#######################################
#######   directory structure  ########
#######################################
installation_path=$(realpath ./)

mkdir -p $installation_path/$pipelineVersion/tools/
mkdir $installation_path/$pipelineVersion/references/
mkdir $installation_path/$pipelineVersion/libraries/
mkdir $installation_path/$pipelineVersion/dbs/

src_path=$installation_path/$pipelineVersion/src
tool_path=$installation_path/$pipelineVersion/tools
dbpath=$installation_path/$pipelineVersion/dbs


#########################################################
#####   Global installation: JAVA  and  python     ######
#########################################################


# java
# optional configuration details at this site:
# https://thishosting.rocks/install-java-ubuntu/  
if type -p java; then
    echo found java executable in PATH
    _java=java
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]];  then
    echo found java executable in JAVA_HOME     
    _java="$JAVA_HOME/bin/java"
else
    echo "no java, installing it just now!"
    sudo apt-get update && apt-get upgrade
    apt-get install default-jdk
fi

if [[ "$_java" ]]; then
    version=$("$_java" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    echo version "$version"
    if [[ "$version" > $min_java ]]; then
        echo version is more than $min_java, OK. 
    else         
        echo version is less than $min_java, please upgrade it 
        exit 1
    fi
fi


## check and install python 2.7 if were necesary
if [ -x "$(command -v python)" ]; then
    echo "python is already instelled in your system"
    # command
else
    echo "Installing python 2.7"
    sudo apt update
    sudo apt upgrade
    sudo apt install python2.7 python-pip
fi



################################################
######   local installations (tools).   ########
################################################

# cromwell
wget $cromwell -P $tool_path 
# wdltool 
wget $wdltool -P $tool_path

# install GATK  (from  docker)
jarname='gatk-package-'$gatkVersion'-local.jar'
gatkzip='gatk-'$gatkVersion'.zip'
folder='gatk-'$gatkVersion
wget $gatklink -P $tool_path/
unzip -p $tool_path/$gatkzip $folder/$jarname > $tool_path/$jarname # extract only the precompiled jarfile


#bwa + reference building
cd $tool_path
wget -O- $bwakit |tar xjf - # get bwa-kit
ln -rs bwa.kit/bwa .            # create link to binary
bwa.kit/run-gen-ref hs37d5      # get reference
bwa.kit/bwa index hs37d5.fa     # build indices
mkdir $installation_path/$pipelineVersion/references/hs37d5/ 
mv hs37d5* $installation_path/$pipelineVersion/references/hs37d5/ # move indices to refernce path

#fastp 
wget $fastp -P $tool_path/
chmod a+x $tool_path/fastp

# dbs (for preprocessing & annotation prouposes)

wget $Mils_1000G_b37_vcf -P $dbpath/       # for preprocessing
wget $Mils_1000G_b37_vcf'.idx' -P $dbpath/

# dbSNP 
wget $dbsnp_b37_vcf -P $dbpath/            # for preprocessing
wget $dbsnp_b37_vcf'.tbi' -P $dbpath/


### For annotation prouposes
$src_path/get_external_dbs.sh $installation_path/$pipelineVersion

