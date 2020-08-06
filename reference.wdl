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
    File bamPath
    File all_regions
    File util
    File count_reads
  }

  command {
    Rscript ${count_reads} ${bamPath} ${all_regions} ""
  }

  runtime {
    time: 10
    cpu: 1
    mem: 10
  }

  output {
     File counts = glob("*.tsv")[0]
  }
}
