workflow TouchWriteHead {
	call touch_write_head
}

task touch_write_head {
	#File GATK
	#File RefFasta
	#File GATK_jar
	#File RefIndex
	#File RefDict
	#File inputBAM
	#File bamIndex
    File localImputFile
	String SampleName
    String contentText

	command {
		tail ${localImputFile} > ${SampleName}.txt
        	}

    runtime {
        docker: "ubuntu"
    }
	output {
	File rawfile = "${SampleName}.txt"
	}

}