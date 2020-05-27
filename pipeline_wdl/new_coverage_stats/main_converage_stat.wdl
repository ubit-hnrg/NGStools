task bam_depth {

File input_bam
File Exon_coords
#String sampleID

###herramientas
String name = basename(input_bam, ".bam")
String toolpath

command <<<

#!/bin/bash
set -e
set -o pipefail

# ok for reduce bam
${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${Exon_coords} > ${name}_exonTSO_reduced.bam

#Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
#que están presentes por una mínima interseccion con el intervalo buscado)
${toolpath}bedtools2/bin/genomeCoverageBed -ibam ${name}_exonTSO_reduced.bam -bga > ${name}_exon.coverage

# Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
${toolpath}bedtools2/bin/intersectBed -a ${Exon_coords} -b ${name}_exon.coverage -loj > ${name}_loj.txt

#esto limita la covertura estrictamente al intervalo paddeado.
# Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
# primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
awk -F"\t" '{print $8"\t"$9"\t"$10"\t"$11"\t"$5"\t"$4"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3}' ${name}_loj.txt > ${name}_loj_sorted_cols.tsv
${toolpath}bedtools2/bin/intersectBed -a ${name}_loj_sorted_cols.tsv -b ${Exon_coords}  > ${name}_loj_exon_filtered.coverage

echo -e 'chr\tstart\tend\tcount\tgeneSymbol\tENSEMBL_ID\texon_number\tstrand\texon_chr\texon_start\texon_end' > header.tsv
cat header.tsv ${name}_loj_exon_filtered.coverage > ${name}_exon_filtered_coverage.tsv
#rm ${name}_loj_exon_filtered.coverage ${name}_loj_sorted_cols.tsv header.tsv ${name}_exon.coverage ${name}_loj.txt

#####coverage statistics cambio la forma del input... 
/home/hnrg/NGStools/pipeline_wdl/qualityControl/coverage_statistics_v1.0.py -i ${name}_exon_filtered_coverage.tsv -g ${name}_global_coverage_statistics.tsv -e ${name}_coverage_statistics_by_exon.tsv -s ${name}

rm ${name}_exonTSO_reduced.bam #${name}_exon_filtered_coverage.tsv
>>>
output {

#### agregar un writelines para generar un archivo .files con el 

File cov_stats_by_exon = "${name}_coverage_statistics_by_exon.tsv"
File glob_cov_stats = "${name}_global_coverage_statistics.tsv"
#String sample_Name = "${name}"

}

}

task make_excel {
String sample_name
File tabla1
String pestana1
#File tabla2
#String pestana2
#File tabla3
#String pestana3

command{
/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${tabla1}:${pestana1} ${sample_name}_exon_coverage_ENS.xlsx
 
}

output {
File reporte_excel = "${sample_name}_exon_coverage_ENS.xlsx"
}
}


task symlink_important_files {
    File output_to_save
    String path_save
    command{
       set -e -o pipefail

       path_out=$(dirname ${path_save})


       cp -L ${output_to_save} $path_out
    }
}
   



workflow quality_control {

String toolpath
File exon_coords
#Array[String] path_save
File list_of_bams


Array[File] bams = read_lines(list_of_bams)
Array[String] path = read_lines(list_of_bams)
Array[Pair[String,File]] bams_x_path = zip(path, bams)

scatter (pairs in bams_x_path)  {

String sample = basename(pairs.left, ".bam")

call bam_depth {
input: 
input_bam = pairs.left,
Exon_coords = exon_coords,
toolpath = toolpath

}

call make_excel {
input:
sample_name = sample,
tabla1 = bam_depth.cov_stats_by_exon,
pestana1 = "exon_coverage_ENS"
#tabla2 = merge_samtools_reports.merged_report, 
#pestana2 = "Alineamiento",
#tabla3 = merge_reports.merged_report,
#pestana3 = "Profundidad-en-libreria"

}

call symlink_important_files {
    input:
    output_to_save = make_excel.reporte_excel,
    path_save = pairs.right
    }


}

#Array[File] bams_stat_depth_global_coverage_stats = bam_depth.glob_cov_stats

#Create a file with a list of the generated bam_depth.glob_cov_stats
#  call CreateFoFN {
#    input:
#      array_of_files = bams_stat_depth_global_coverage_stats,
#      fofn_name = Tso_name
     
#  }

Array[File] reportes_salidas = ["${make_excel.reporte_excel}"]





}