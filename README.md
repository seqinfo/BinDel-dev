# BinDel

Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
All rights reserved, unauthorised usage and distribution are prohibited.
Contact: priit.palta@gmail.com

#### Abstract

The scientific software focuses on detecting rare (occurring in five or fewer people in 10,000) clinically relevant pathogenic microdeletions from low-coverage NIPT WGS data. 
However, the software is not limited to microdeletion detection and is developed with the idea of detecting any difference from the reference set. Detection possibility includes the detection of full chromosome aneuploidies or monosomies. When a female fetus reference set is used, then the software supports 45,X detection. Moreover, BinDel can focus on any predefined subregion in the chromosome with varying bin lengths predefined in the `coordinates.bed` file.

#### Algorithm
The algorithm applies for each bin bin-based [GC% correct](https://dx.doi.org/10.1038%2Fs41598-017-02031-5), normalises bins by bin length and a sample total read count. Next, the software applies [PCA-based normalisation](https://doi.org/10.1038/gim.2018.32), calculates Z-scores based on the reference bins, applies Z-score normalisation based on the sample region mean read count, calculates Mahalanobis distance from the reference and converts them via Chi-Square distribution to high-risk probabilities.


# Manual

## Alignment
BinDel requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Coordinates
BinDel requires a `.bed` file with predefined coordinates. The coordinate file describes:
* Chromosomes (column `chr`)
* *bins* and bin lengths (bins can have varying length) are defined by columns `start` and `end`.
* Names of the subregions of interest to analyse (column `focus`). It can be a whole chromosome or a subregion.

**Note 1:** Columns `chr`, `start` and `end` must uniquely define each region, e.g. `.bed` file must not contain duplicates. Column `focus` is the name of the region of interest, which means that this column is used for grouping bins. **Having duplicates in .bed leads to anomalies in final high-risk probabilities**.

**Note 2:** GC% correct depends on the number of regions of interest. E.g. if only, for example, chromosome 2 is in the analysis, it can affect the risk scoring compared to having all chromosomes in the analysis.

**Note 3:** **For the detection of microdeletions it is important to prefilter regions with high variance.** For example, we observed that Prader-Willi (PWS) and Angelman syndrome (AS) related region has high variance in the beginning of the region. In the development of the tool, only coordinates `24500001 - 27800000` were used for AS/PWS detection.

Example of the `.bed`:

| chr  | start | end | focus |
| ------------- | ------------- | ------------- | ------------- |
| chr1  | 1  | 100000 | 1p36 |
| chr1  | 100001  | 200000 | 1p36 |
| chr1  | 200001  | 300000 | 1p36 |
| chr1  | 300001  | 400000 | 1p36 |
| ...  | ...  | ... | ... |
| ...  | ...  | ... | ... |
| chr22  | 17200001  | 17300000 | chr22 |
| chr22  | 19100001  | 19200000 | DiGeorge |
| chr22  | 19200001  | 19300000 | DiGeorge |
| ...  | ...  | ... | ... |
| ...  | ...  | ... | ... |
| chr22  | 50600001  | 50700000 | chr22 |
| chr22  | 50700001  | 50800000 | chr22 |

#### Automated creation of the `.bed`:

Given a file `locations.info.tsv` (column `length` sets bin size, e.g in here bin size is 100k):
```TSV
chr	start	end	focus	length
...
chr5	1	48800000	Cri-du-chat	100000
chr5	48800001	181538259	chr5	100000
chr6	1	170805979	chr6	100000
...

```
Run the following Python script:
```
python3 dividebins.py --infile locations.info.tsv --outfile coordinates.bed
```
The script creates the file `coordinates.bed`, which can be used in the reference file creation. Note that this example excludes chromosome X and Y. However, there is no software side limitation of including chromosomes X and Y in the analysis. Allosomes are processed in the same way as autosomes, and therefore if 45,X is the subject of interest, only female fetus reference group should be used, and in the analysable sample, only a female fetus sample should be used. 
## Installation
```R
# In R:
install.packages("devtools")
devtools::install_github("seqinfo/BinDel")
```
## Reference creation
BinDel requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires the existence of `coordinates.bed` file, which defines the subregions in the genome to analyse.
```R
# In R:
bindel::create_reference("path/folder/bams", "path/coordinates.bed", "name_of_the_output")
```
## Usage
```R
# In R:
bindel::infer_normality("path/bam.bam", "path/reference.tsv")
```
Note, if the reference file has fewer samples than the default number of PCA components to be used in the normalisation, set the parameter `nComp <= number of reference samples` or turn off PCA normalisation.

```R
# In R:
bindel::infer_normality("path/bam.bam", "path/reference.tsv", nComp = less_than_n_samples_in_reference)
```

For more information about possible parameters, please consult the function documentation:
```R
# In R:
?bindel::infer_normality
```

## Output
`bindel::infer_normality("path/bam.bam", "reference_location.tsv")` outputs by default two scientific files:
1. `.png` containing normalised Z-score for each bin in each subregion. These Z-scores are the basis of the high-risk probability calculation. These figures also illustrate the reference group Z-scores.
2. `.tsv` summary file for each subregion.
