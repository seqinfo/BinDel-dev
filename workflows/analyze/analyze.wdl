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
    Array[File] scores = analyze.results 
    Array[File] plots = analyze.plots 
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
     File results = glob("*.results.tsv")[0]
     File plots = glob("*.pdf")[0]     
  }
}
