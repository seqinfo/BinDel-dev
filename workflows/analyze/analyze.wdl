version 1.0

workflow Main {

    input {
        Array[String] bam_locations
        File coordinates
    }

    scatter (bamPath in bam_locations) {
        call analyze {
            input:
                bamPath = bamPath
        }
    }

    output {
        Array[File] scores = analyze.result
        Array[File] plots = analyze.plot
    }
}

task analyze {

    input {
        File bamPath
        File reference
    }

    command {
        set -e
        Rscript $analyser ${bamPath} ${reference}
    }

    runtime {
        time: 50
        cpu: 1
        mem: 50
    }
    output {
        File result =  basename(bamPath) + ".tsv"
        File plot = basename(bamPath) + ".png"
    }
}
