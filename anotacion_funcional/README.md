## Pipeline para procesar raw data de NGS. 
* Pipeline of 14 steps of funtional annotation. 

### Notes: 
* It´s strongly recommended to use the provided cromwell's configuration file in order to enable softlink usage. 
* Tested in Cromwell 33-1 version. 

## Installation steps. 

# Requerimentos: 
* VCF file.
* wdl configuration file (json format). This file contains paths to input data as well as database and reference genomes location paths.  

## usage examples.
* from solari server run: 

```bash 
sudo java -Dconfig.file=./archi.conf -jar /tools/cromwell/cromwell-33.1.jar run anotaciones_actualizadas.wdl -i inputs_anotaciones_act.json 
```
