# BINDEL

The software focuses on detecting rare (occurring in five or fewer people in 10,000) clinically relevant pathogenic microdeletions from low-coverage NIPT WGS data. 
However, the software is not limited to microdeletion detection and is developed with the idea of detecting any difference from the reference set. Detection possibility includes the detection of full chromosome aneuploidies or monosomies. When a female fetus reference set is used, then the software supports 45,X detection.

For examples and usage, please consult with the Cromwell workflow in the corresponding GitHub [repository](https://github.com/seqinfo/PPDxWorkflow).

Software is developed for GRCh38 genomic version.

# Manual

## Alignment
BINDEL requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Coordinates
BINDEL requires a `.bed` file with predefined coordinates. Coordinate file describes:
* Chromosomes (column `chr`)
* *bins* and bin lengths (bins can have varying length) which are defined by columns `start` and `end`.
* Names of the subregions of interest to analyse (column `focus`). Can be a whole chromosome or a subregion.

**Note:** Columns `chr`, `start` and `end` must uniquely define each region, e.g `.bed` file must not contain duplicates. Column `focus` is the name of the region of interest, which means that this column is used for grouping bins. **Having duplicates in .bed leads to anomalies in final high-risk probabilites**.

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

## Installation
```R
# In R:
install.packages("devtools")
devtools::install_github("seqinfo/BINDEL/")
```
## Reference creation
```R
# In R:
bindel::create_reference("path/folder/bams", "path/coordinate_bed.bed", "name_of_the_output")
```
