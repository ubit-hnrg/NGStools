#!/bin/bash
for vcf in "$@"
do
sed -i -e 's|Veredict=Uncertain significance|Veredict=Uncertain_significance|g' $vcf
sed -i -e 's|Veredict=Likely pathogenic|Veredict=Likely_pathogenic|g' $vcf
sed -i -e 's|Veredict=Likely benign|Veredict=Likely_benign|g' $vcf
done