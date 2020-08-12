version 1.0

workflow Main {

  input {
    Array[String] bam_locations
    File coordinates
    File coordinates_of_interest
  }

  scatter (bamPath in bam_locations) {
    call analyze {
      input:
        bamPath = bamPath,
        coordinates = coordinates
    }

    call visualize {
      input:
        tsv = analyze.results,
        coordinates_of_interest = coordinates_of_interest
    }
  }

  output {
    File scores = analyze.results 
    File plots = visualize.plots 
  }
}

task analyze {

  input {
    File bamPath
    File coordinates
    File reference
    }

  command {
    set -e
    Rscript $analyser ${bamPath} ${coordinates} ${reference}
  }

  runtime {
    time: 10
    cpu: 1
    mem: 10
  }
    output {
     File results = glob("*.tsv")[0]
  }
}

task visualize {

  input {
    File tsv
    File coordinates_of_interest
    }

  command {
    Rscript $visualizer ${tsv} {coordinates_of_interest}
  }

  runtime {
    time: 10
    cpu: 1
    mem: 10
  }

  output {
     File plots = glob("*.pdf")[0]
  }
}
