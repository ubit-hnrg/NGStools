#!/bin/bash
###########################################
############# Sacar herramientas de diferentes docker
############# BWA Y GATK
############# SE USA DOCKER "GENOMES ON THE CLOUD" y GATK4 PARA USAR HERRAMIENTAS EMPAQUETADAS
#################################################

tool_path=$(realpath ~/tools/)
mkdir $tool_path
gatk_symlink='gatk.jar'

#### ID genomes on the CLOUD
#name=broadinstitute/genomes-in-the-cloud
#tag=2.3.1-1512499786
#dockerid1=$(docker inspect --format="{{.Id}}" $name:$tag |cut -f2 -d":"|cut -c1-12)
#sudo docker pull $name:$tag



###llamar a docker genomes on the cluod y sacar BWA
#docker run -v $tool_path:/herramientasdocker $dockerid1 cp -r bgzip /herramientasdocker 


#### ID GATK4
name='broadinstitute/gatk'
tag='4.0.6.0'
dockerid2=$(docker inspect --format="{{.Id}}" $name:$tag |cut -f2 -d":"|cut -c1-12)
echo $dockerid2
sudo docker pull $name:$tag

### llama GATK$ y saca el empaquetado de herramientas gatk-"version"-local.jar

# parche para dereferenciar el symlink gatk.jar
jarname=$(docker run $dockerid2 ls -lrth $gatk_symlink |cut -f2 -d">"|cut -f2 -d' '|sed -e "s/\r//g") 

#copy gatk-package-[VERSION]-local.jar
docker run -v $tool_path:/dockertools/ $dockerid2 cp $jarname /dockertools/
sudo chown -R $USER:$USER $tool_path
