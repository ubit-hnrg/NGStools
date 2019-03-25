workflow test {
    File in
    call create_dummyfile { input: in=in }

    call symlink_a_lo_guapo {
        input:
        output_to_save = create_dummyfile.dummy_outfile
    }
}

task create_dummyfile {
    File in
    command {
        cat ${in} > blabla.txt
    }
    
    output{
        File dummy_outfile = "blabla.txt"
    }
}

task mkdir_samplename {
    String samplename
    command{
        mkdir ${samplename} /home/hnrg/resultsHNRG/
    }
}


task symlink_a_lo_guapo {
    File output_to_save
    String path = /home/ariel/Documents/
    command{
        ln -s ${output_to_save} ${path}
    }
}
