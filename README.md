# SyndromeDetector


(Genomic coordinates from A534N .bed.)


Usage:

```mkdir "temp"```

```Rscript count_reads.R "60569250_S11.bam" "coordinates/chr15.bed" "temp"```

```Rscript create_ref.R "temp" "reference.tsv" ```

```Rscript analyse_bam.R "60569250_S11.bam" "coordinates/chr15.bed" "reference.tsv" ```

```Rscript visualize.R "60569250_S11.results.tsv" "60569250_S11" "coordinates/as_pws.bed" "AS/PWS" ```

