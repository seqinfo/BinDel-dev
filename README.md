# BinDel: software tool for detecting clinically significant microdeletions in low-coverage WGS-based NIPT samples
BinDel is distributed under the Attribution-NonCommercial-ShareAlike 4.0 International ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)) license.

Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
All rights reserved, unauthorised usage and distribution are prohibited.
Contact: priitpaluoja@gmail.com/priit.palta@gmail.com



Here we present the BinDel, a novel region-aware microdeletion detection software package developed to infer clinically relevant microdeletion risk in low-coverage whole-genome sequencing NIPT data. 

Our [paper](https://doi.org/10.1101/2022.09.20.22280152) describes the BinDel algorithm and how it was tested. We quantified the impact of sequencing coverage, fetal DNA fraction, and region length on microdeletion risk detection accuracy. We also estimated BinDel accuracy on known microdeletion samples and clinically validated aneuploidy samples. 


## Installation

```R
# In R:
install.packages("devtools") # Skip this line if devtools is already installed
devtools::install_github("seqinfo/BinDel")
```
<details><summary>Note in case of XML dependency mismatch</summary>
<p>

 We encountered an error of unable to install R package due to XML dependency mismatch on one of the test computers. We solved it by installing `r-xml` with `conda install r-xml`). Our solution was based on [this](https://stackoverflow.com/questions/37035088/unable-to-install-r-package-due-to-xml-dependency-mismatch).

```
Using libxml2.*
checking for gzopen in -lz... yes
checking for xmlParseFile in -lxml2... no
checking for xmlParseFile in -lxml... no
configure: error: "libxml not found"
ERROR: configuration failed for package 'XML'
 ```
 ```
 Error in `(function (command = NULL, args = character(), error_on_status = TRUE, ...`:
! System command 'R' failed
---
Exit status: 1
stdout & stderr: <printed>
---
Type .Last.error to see the more details.
Warning messages:
1: In i.p(...) : installation of package 'XML' had non-zero exit status
2: In i.p(...) :
  installation of package 'restfulr' had non-zero exit status
3: In i.p(...) :
  installation of package 'rtracklayer' had non-zero exit status
4: In i.p(...) :
  installation of package 'BSgenome' had non-zero exit status
5: In i.p(...) :
  installation of package 'BSgenome.Hsapiens.UCSC.hg38' had non-zero exit status
```
 
</p>
</details>


## Alignment/Mapping
BinDel requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Usage
### Reference creation
BinDel requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires known euploid NIPT samples in `.bam` format and
a [file](example/bins.bed) defining each genomic bin to use in the analysis.

```R
# In R:
bindel::create_reference("path/folder/bams", "example/bins.bed", "name_of_the_reference_file")
```

<details><summary> Click here to see how to define genomic bins.</summary>
<p>

Given a file [`example/locations.info.tsv`](example/locations.info.tsv) describing bin lengths (column `length`) for each region of interest, the following Python [script](dividebins.py) bins the input file to [`example/bins.bed`](example/bins.bed):
 
```
python3 dividebins.py --infile example/locations.info.tsv --outfile example/bins.bed
```
The script creates the file [`example/bins.bed`](example/bins.bed), which can be used in the reference file creation.

<details><summary>Notes</summary>
<p>

**Note 1:** Columns `chr`, `start` and `end` must uniquely define each region, e.g. `.bed` file must not contain duplicates. Column `focus` is the name of the region of interest, which means that this column is used for grouping bins. **Having duplicates in .bed leads to anomalies in final high-risk probabilities**.

**Note 2:** GC% correct depends on the number of regions of interest. E.g. if only, for example, chromosome 2 is in the analysis, it can affect the risk scoring compared to having all chromosomes in the analysis.
</p>
</details>

</p>
</details>


### Running BinDel
```R
# In R:
bindel::infer_normality("sample.bam", "name_of_the_reference_file")
```

<details><summary>Note</summary>
<p>

If the reference file has fewer samples than the default number of PCA components to be used in the normalisation, set the parameter `nComp` to lower than number of reference samples used.

```R
# In R:
bindel::infer_normality("path/bam.bam", "path/reference.tsv", nComp = less_than_n_samples_in_reference)
```
</p>
</details>
