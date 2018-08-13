workflow hpc {
	call HPCaller
}

task HPCaller {
	File RefFasta
	File RefIndex
	File RefDict

	File inputBAM
	File bamIndex

	String SampleName

	command {
		gatk --java-options "-Xmx4g" HaplotypeCaller \
		-R ${RefFasta} \
		-I ${inputBAM} \
		-O ${SampleName}.raw.indels.snps.g.vcf  \
		-L 20 \
        -new-qual\
        -ERC GVCF
	}
	runtime {
#		docker: "broadinstitute/gatk"
		docker: "1171e1d71355"

	}

	output {
	File rawVCF = "${SampleName}.raw.indels.snps.g.vcf"
	}
}


