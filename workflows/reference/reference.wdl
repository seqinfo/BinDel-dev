version 1.0

workflow Main {
  input {
    Array[String] bam_locations
  }

  scatter (bamPath in bam_locations) {
    call bin {
      input:
        bamPath = bamPath
    }
  }

  call merge {
    input:
      counts = bin.counts
  }

  output {
    File reference = merge.reference 
  }
}


task bin {

  input {
    File bamPath
    File regions
    }

  command {
    Rscript ../R/count_reads.R ${bamPath} ${regions} ""
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


task merge {

  input {
    File counts
    }

  command {
    Rscript ../R/create_ref.R ${counts} "reference.tsv"
  }

  runtime {
    time: 20
    cpu: 1
    mem: 10
  }

  output {
     File reference = "reference.tsv"
  }
}
