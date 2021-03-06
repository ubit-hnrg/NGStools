workflow main_workflow {

File GatheredBamFile
File TSO_padded_bedfile
File exon_coord_inTSO_bedfile

call coverBed{
    input:
    GatheredBamFile=GatheredBamFile
    TSO_padded_bedfile=TSO_padded_bedfile    
}

task coverBed{
    File GatheredBamFile
    File TSO_padded_bedfile
    String toolpath
    String coverage_by_position_suffix
    command{
        ${toolpath}/bedtools2/bin/coverageBed -a %s ${TSO_padded_bedfile} ${GatheredBamFile} -d > ${coverage_by_position_suffix}.txt
    }
}

toolpath=/home/bitgenia/samples/HNRG-pipeline-V0.1/tools
$toolpath/bedtools2/bin/coverageBed -a %s ${TSO_padded_bedfile} ${GatheredBamFile} -d > ${coverage_by_position_suffix}.txt


sampleID=$1
BAM='/home/hnrg/resultsHNRG/'$sampleID'/'$sampleID'.bam'
EXON_coords=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_in_TSO_V1_padded.bed
# ok for reduce bam
$toolpath/bedtools2/bin/intersectBed -a $BAM -b $EXON_coords > $sampleID'_exonTSO_reduced.bam'

#Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
#que están presentes por una mínima interseccion con el intervalo buscado)
$toolpath/bedtools2/bin/genomeCoverageBed -ibam $sampleID'_exonTSO_reduced.bam' -bga > $sampleID'_exon.coverage'

# Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
$toolpath/bedtools2/bin/intersectBed -a EXON_coords -b $sampleID'_exon.coverage' -loj > loj.txt


#esto limita la covertura estrictamente al intervalo paddeado.
# Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
# primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
awk -F"\t" '{print $1"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$2"\t"$3}' loj.txt > loj_sorted_cols.tsv
$toolpath/bedtools2/bin/intersectBed -a loj_sorted_cols.tsv -b EXON_coords  > $sampleID'_loj_exon_filtered.coverage'

echo -e 'chr\tstart\tend\tcount\tgeneSymbol\texon_number\texon_start\texon_end' > header
cat header $sampleID'_loj_exon_filtered.coverage' > $sampleID'_exon_filtered_coverage.tsv'
rm $sampleID'_loj_exon_filtered.coverage loj.txt' loj_sorted_cols.tsv header $sampleID'_exon.coverage'


