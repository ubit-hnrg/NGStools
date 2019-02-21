#!/bin/bash
# web resources : 
#tools and references
cromwell='https://github.com/broadinstitute/cromwell/releases/download/37/cromwell-37.jar'
wdltool='https://github.com/broadinstitute/wdltool/releases/download/0.14/wdltool-0.14.jar'
bwakit='http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download'
fastp='http://opengene.org/fastp/fastp'

#libraries 
wgs_interval_list='http://???'
wes_interval_list='http://???'

# dbs 
Mils_1000G_b37_vcf='http://???'         #Mills&1000G
Mils_1000G_b37_vcf_idx='http://???'
dbsnp_b37_vcf='http://???'              #dbSNP
dbsnp_b37_vcf_idx='http://???' 





# parameters
docker_gatk_name='broadinstitute/gatk'
tag='4.0.6.0'



########################################
#####    Global installations. #########
########################################

installation_path=$(realpath ./)

# java
# optional configuration details at this site:
# https://thishosting.rocks/install-java-ubuntu/  
sudo apt-get update && apt-get upgrade
apt-get install default-jdk

## install docker

# check python 2.7


########################################
######   local installations.   ########
########################################

# directory structure
mkdir -p $installation_path/HNRG-pipeline/tools/
mkdir $installation_path/HNRG-pipeline/references/
mkdir $installation_path/HNRG-pipeline/libraries/
mkdir $installation_path/HNRG-pipeline/dbs/

# cromwell
cd $installation_path/HNRG-pipeline/tools/
wget $cromwell

# wdltool 
cd $installation_path/HNRG-pipeline/tools/
wget $wdltool


# install GATK  (from  docker)
#gatk_symlink='gatk.jar'
tool_path=$installation_path/HNRG-pipeline/tools/
dockerid2=$(docker inspect --format="{{.Id}}" $docker_gatk_name:$tag |cut -f2 -d":"|cut -c1-12)
echo 'downloading' $dockerid2
sudo docker pull $docker_gatk_name:$tag
#jarname=$(docker run $dockerid2 ls -lrth $gatk_symlink |cut -f2 -d">"|cut -f2 -d' '|sed -e "s/\r//g") # parche para dereferenciar el symlink gatk.jar
jarname='gatk-package-'$tag'-local.jar'
docker run -v $tool_path:/dockertools/ $dockerid2 cp $jarname /dockertools/ #copy gatk-package-[VERSION]-local.jar
sudo chown -R $USER:$USER $tool_path  #change permisions. 


#bwa + reference building
cd $installation_path/HNRG-pipeline/tools/
wget -O- $bwakit |tar xjf - # get bwa-kit
ln -rs bwa.kit/bwa .            # create link to binary
bwa.kit/run-gen-ref hs37d5      # get reference
bwa.kit/bwa index hs37d5.fa     # build indices
mkdir $installation_path/HNRG-pipeline/references/hs37d5/ 
mv hs37d5* $installation_path/HNRG-pipeline/references/hs37d5/ # move indices to refernce path


#fastp 
cd $installation_path/HNRG-pipeline/tools/
wget $fastp -P $installation_path/HNRG-pipeline/tools/
chmod a+x $installation_path/HNRG-pipeline/tools/fastp

# libraries /references
# download interval lists and bedfiles for wes and wgs intervals.
wget $wgs_interval_list -P $installation_path/HNRG-pipeline/libraries/
wget $wes_interval_list -P $installation_path/HNRG-pipeline/libraries/


# dbs (for preprocessing & annotation prouposes)

wget $Mils_1000G_b37_vcf -P $installation_path/HNRG-pipeline/dbs/       # for preprocessing
wget $Mils_1000G_b37_vcf_idx -P $installation_path/HNRG-pipeline/dbs/

wget $dbsnp_b37_vcf -P $installation_path/HNRG-pipeline/dbs/            # for preprocessing
wget $dbsnp_b37_vcf_idx -P $installation_path/HNRG-pipeline/dbs/
