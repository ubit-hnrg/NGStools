
# anotar primero con anovar. 


#cambio luego el nombre de las samples Numericas (sino SnpSift lo interpreta como entero)
cp TSO20192803_only_kit_annovar.hg19_multianno.vcf TSO20192803_only_kit_renamed_annovar.hg19_multianno.vcf
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 ;do echo $i;sed -i -e"s/$i/ID$i/g" TSO20192803_only_kit_renamed_annovar.hg19_multianno.vcf;done

# #filtro con SnpSift pidiendo para cada muestra que un vcf que no contenga 0/0 ni ./. de ella (si de las otras) 
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;do echo $i; cat TSO20192803_only_kit_renamed_annovar.hg19_multianno.vcf | java -jar ~/HNRG-pipeline-V0.1/tools/SnpSift.jar filter "GEN[ID$i].GT != './.' && GEN[ID$i].GT != '0/0'" > /home/hnrg/resultsHNRG/$i/$i'_only_kit_renamed_annovar.hg19_multianno.vcf';done

# Ahora puedo correr el reporte de vcf en estos archivos tal cual está en el pipeline. 
#---> Pendiente


# Ahora tengo que filtrar la tabla con estas mismas variantes. 


# y tambien tengo que filtrar esa tabla con un subset de genes.


#y con eso estaríamos. !