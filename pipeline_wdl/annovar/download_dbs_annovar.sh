#refGene
##wget http://www.openbioinformatics.org/annovar/download/hg19_refGene.txt.gz
##wget http://www.openbioinformatics.org/annovar/download/hg19_refGeneMrna.fa.gz
##wget http://www.openbioinformatics.org/annovar/download/hg19_refGeneVersion.txt.gz

#esp6500siv2_all
##wget http://www.openbioinformatics.org/annovar/download/hg19_esp6500siv2_all.txt.gz
##wget http://www.openbioinformatics.org/annovar/download/hg19_esp6500siv2_all.txt.idx.gz

#1000g2015aug
##wget http://www.openbioinformatics.org/annovar/download/hg19_1000g2015aug.zip

#avsnp150
#wget http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.idx.gz

#dbnsfp35a
#wget http://www.openbioinformatics.org/annovar/download/hg19_dbnsfp35a.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_dbnsfp35a.txt.idx.gz

#clinvar_20180603
#wget http://www.openbioinformatics.org/annovar/download/hg19_clinvar_20180603.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_clinvar_20180603.txt.idx.gz

#gnomad_genome
#wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.idx.gz

#dbscsnv11
#wget http://www.openbioinformatics.org/annovar/download/hg19_dbscsnv11.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_dbscsnv11.txt.idx.gz
                                                                         -

#rmsk
##wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz

#ensGene
#wget http://www.openbioinformatics.org/annovar/download/hg19_ensGene.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_ensGeneMrna.fa.gz

#knownGene
##wget http://www.openbioinformatics.org/annovar/download/hg19_knownGene.txt.gz
##wget http://www.openbioinformatics.org/annovar/download/hg19_kgXref.txt.gz
##wget http://www.openbioinformatics.org/annovar/download/hg19_knownGeneMrna.fa.gz


# intervar_20180118
#wget http://www.openbioinformatics.org/annovar/download/hg19_intervar_20180118.txt.gz 
#wget http://www.openbioinformatics.org/annovar/download/hg19_intervar_20180118.txt.idx.gz

# exac03
#wget http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.idx.gz

# genomad_exome
#wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_exome.txt.gz
#wget http://www.openbioinformatics.org/annovar/download/hg19_gnomad_exome.txt.idx.gz

annovar_urlpath='http://www.openbioinformatics.org/annovar/download/hg19_'
golden_urlpath='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database'

#for db in wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75 wgEncodeBroadChipSeqPeaksGm12878H3k4me1 wgEncodeRegDnaseClustered  # all these dbs are only in hg18 version
for db in tfbsConsSites cytoBand wgRna targetScanS genomicSuperDups dgvMerged gwasCatalog wgEncodeRegTfbsClustered
do
    wget $golden_urlpath/$db.txt.gz
    #wget $golden_urlpath$db.txt.idx.gz
    gzip -d $db.txt.gz
    #gzip -d $db.txt.idx.gz
    mv $db.txt hg19_$db.txt
    #mv $db.txt.idx.gz hg19_$db.txt.idx.gz
done


