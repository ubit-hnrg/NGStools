workflow test {
    File in
    call cat { input: in=in }
}

task cat {
    File in
    command {
        cat ${in}
    }
}
