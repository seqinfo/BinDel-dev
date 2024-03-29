# BinDel: software tool for detecting clinically significant microdeletions in low-coverage WGS-based NIPT samples
BinDel is distributed under the Attribution-NonCommercial-ShareAlike 4.0 International ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)) license.

Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
All rights reserved, unauthorised usage and distribution are prohibited.
Contact: priitpaluoja@gmail.com/priit.palta@gmail.com



Here we present the BinDel, a novel region-aware microdeletion detection software package developed to infer clinically relevant microdeletion risk in low-coverage whole-genome sequencing NIPT data. 

Our [paper](https://doi.org/10.1101/2022.09.20.22280152) describes the BinDel algorithm and how it was tested. We quantified the impact of sequencing coverage, fetal DNA fraction, and region length on microdeletion risk detection accuracy. We also estimated BinDel accuracy on known microdeletion samples and clinically validated aneuploidy samples. 


## Installation
#### Docker: [![Docker Repository on Quay](https://quay.io/repository/priitpaluoja/bindel/status "Docker Repository on Quay")](https://quay.io/repository/priitpaluoja/bindel)
<details><summary>Installation on Ubuntu 22.04</summary>
<p>
 
The following is tested with [ubuntu-22.04.1-live-server-amd64](https://releases.ubuntu.com/22.04/).

#### Install R as shown in [DigitalOcean](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04). [From DigitalOcean](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-22-04):
```bash
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
sudo apt update
sudo apt install --no-install-recommends r-base
```
#### Install BinDel dependencies and [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
```bash
sudo apt -y install r-cran-devtools r-bioc-biostrings r-cran-dplyr r-bioc-genomicalignments r-bioc-genomicranges r-cran-ggplot2  r-bioc-iranges r-cran-magrittr r-cran-purrr r-cran-readr r-bioc-rsamtools r-cran-stringr  r-cran-tidyr git r-bioc-bsgenome  libcairo2-dev libxt-dev
```
#### Install BSgenome.Hsapiens.UCSC.hg38 and BinDel
```R
sudo -i R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
devtools::install_github("seqinfo/BinDel", upgrade = "never")
```
</p>
</details>


<details><summary>Installation on Windows 10</summary>
<p>

1. Install [R](https://cran.r-project.org/bin/windows/base/).
2. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
3. Install devtools and BinDel in R:
```R
# In R:
install.packages("devtools") 
devtools::install_github("seqinfo/BinDel")
```
</p>
</details>


## Alignment/Mapping
BinDel requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Usage
### Reference creation
BinDel requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires known euploid NIPT samples in `.bam` format and
a [coordinates file](example/locations.info.tsv) describing bin lengths and region of interests.

```R
# In R:
BinDel::write_reference(c("sample1.bam", "sample2.bam"), "coordinates.tsv", "reference.gz")
```

### Running BinDel
```R
# In R:
BinDel::infer_normality("sample.bam", "reference.gz")
```
