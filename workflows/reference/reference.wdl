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
    set -e
    Rscript $count_reads ${bamPath} ${regions} ""
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
    Array[File] counts
    File file1 = counts[0]
    }

  command {
    head -1 ${file1} > reference.tsv; tail -n+2 -q ${counts} >> reference.tsv
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
