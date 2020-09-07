version 1.0

workflow Main {

  input {
    Array[String] bam_locations
    File coordinates   
  }

  scatter (bamPath in bam_locations) {
    call analyze {
      input:
        bamPath = bamPath,
        coordinates = coordinates
    }

    call visualize {
      input:
        tsv_results = analyze.results,
        tsv_segments = analyze.segments,
        coordinates_of_interest = coordinates
    }
  }

  output {
    Array[File] scores = analyze.results 
    Array[File] segments = analyze.segments
    Array[File] plots = visualize.plots 
    Array[File] genes = visualize.genes 
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
     File results = glob("*.results.tsv")[0]
     File segments = glob("*.segments.tsv")[0]
  }
}

task visualize {

  input {
    File tsv_results
    File tsv_segments
    File coordinates_of_interest
    }

  command {
    Rscript $visualizer ${tsv_results} ${tsv_segments} ${coordinates_of_interest}
  }

  runtime {
    time: 10
    cpu: 1
    mem: 10
  }

  output {
     File plots = glob("*.pdf")[0]
     File genes = glob("*.genes.tsv")[0]
  }
}
