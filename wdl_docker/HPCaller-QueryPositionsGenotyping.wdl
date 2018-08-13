workflow hpc {
    File RefFasta
	File RefIndex
	File RefDict
    String SampleName
    File queryList

	call HPCaller {
        input:
        prefix = SampleName,
        RefFasta = RefFasta,
        RefDict = RefDict,
        RefIndex = RefIndex,
        queryList = queryList
        
    }

    call genotypeGVCF{
        input:
        prefix = SampleName,
        RefFasta = RefFasta,
        RefDict = RefDict,
        RefIndex = RefIndex,
        queryList = queryList,
        rawGvcf = HPCaller.rawGvcf,
        rawGvcfIndex = HPCaller.rawGvcfIndex

    }

}

task HPCaller {   
    File RefFasta
	File RefIndex
	File RefDict
    File queryList
	File inputBAM
	File bamIndex
    String prefix

	command {
		gatk --java-options "-Xmx4g" HaplotypeCaller \
		-R ${RefFasta} \
		-I ${inputBAM} \
		-O ${prefix}.raw.indels.snps.g.vcf  \
		-L ${queryList} \
        -new-qual\
        -ERC BP_RESOLUTION
	}
	runtime {
#		docker: "broadinstitute/gatk"
		docker: "1171e1d71355"  # #este id es el de la imagen del docker. Ponerlo explicitamente y no con el nombre de la imagen, evita que, si
                                # GATK hizo un update respecto al que vos tenes localmente instalado, te baje y genere una copia del actual. 

	}

	output {
	    File rawGvcf = "${prefix}.raw.indels.snps.g.vcf"
        File rawGvcfIndex = "${prefix}.raw.indels.snps.g.vcf.idx"
	}

}

task genotypeGVCF{
    File RefFasta
	File RefIndex
	File RefDict
    File queryList
    File rawGvcf
    File rawGvcfIndex
    String prefix

    command{
        gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${RefFasta} \
        -V ${rawGvcf} \
        -L ${queryList} \
        -O ${prefix}.raw.vcf
    }
	runtime {
#		docker: "broadinstitute/gatk"
		docker: "1171e1d71355" 
	}

	output {
	    File rawVCF = "${prefix}.raw.vcf"
	}

}




