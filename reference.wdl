version 1.0

workflow Main {
  input {
    Array[String] bam_locations
  }

  scatter (bamPath in bam_locations) {
    call count {
      input:
        bamPath = bamPath
    }
  }

  output {
    Array[File] counts = count.counts 
  }
}

task count {
  input {
    String bamPath
    String all_regions
    String count_reads
    String util
  }

  command {
    Rscript ${count_reads} ${bamPath} ${all_regions} ""
  }

  runtime {
    time: 10
    cpu: 1
    memory: 10
  }

  output {
     File counts = glob("*.tsv")[0]
  }
}
