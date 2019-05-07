#!/bin/bash
# convert multianno into annovar vcf db (sorted)
input=/home/hnrg/resultsHNRG/1711242/1711242_TSO_renamed_one_sample.hg19_multianno.txt
out=17111242_intervarDB.vcf
python NGStools/pipeline_wdl/anotacion_funcional/create_InterVarDB.py -i=$input -o $out

# index
bgzip $out
tabix $out.gz

