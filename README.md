# BINDEL
#### Abstract

The scientific software focuses on detecting rare (occurring in five or fewer people in 10,000) clinically relevant pathogenic microdeletions from low-coverage NIPT WGS data. 
However, the software is not limited to microdeletion detection and is developed with the idea of detecting any difference from the reference set. Detection possibility includes the detection of full chromosome aneuploidies or monosomies. When a female fetus reference set is used, then the software supports 45,X detection. Moreover, BINDEL can focus on any predefined subregion in the chromosome with varying bin lengths predefined in the `coordinates.bed` file.

#### Algorithm
The algorithm applies for each bin bin-based [GC% correct](https://dx.doi.org/10.1038%2Fs41598-017-02031-5), normalises bins by bin length and a sample total read count. Next, the software applies [PCA-based normalisation](https://doi.org/10.1038/gim.2018.32), calculates Z-scores based on the reference bins, applies Z-score normalisation based on the sample region mean read count, calculates Mahalanobis distance from the reference and converts them via Chi-Square distribution to high-risk probabilities.


# Manual

## Alignment
BINDEL requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Coordinates
BINDEL requires a `.bed` file with predefined coordinates. The coordinate file describes:
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

Given a file `locations.info.tsv`:
```TSV
chr	start	end	focus	length
chr1	10001	12780116	1p36 (1/3)	100000
chr1	12780117	23600000	1p36 (2/3)	100000
chr1	23600001	27600000	1p36 (3/3)	100000
chr1	27600001	248956422	chr1	100000
chr2	1	196060396	chr2	100000
chr2	196060397	204342216	2q33.1	100000
chr2	204342217	242193529	chr2	100000
chr3	1	192600000	chr3	100000
chr3	192600001	198295559	3q29	100000
chr4	1	4500000	Wolf-Hirschhorn	100000
chr4	4500001	190214555	chr4	100000
chr5	1	48800000	Cri-du-chat	100000
chr5	48800001	181538259	chr5	100000
chr6	1	170805979	chr6	100000
chr7	1	72700000	chr7	100000
chr7	72700001	77900000	Williams-Beuren	100000
chr7	77900001	159345973	chr7	100000
chr8	1	116700000	chr8	100000
chr8	116700001	126300000	Langer-Giedion	100000
chr8	126300001	145138636	chr8	100000
chr9	1	138394717	chr9	100000
chr10	1	133797422	chr10	100000
chr11	1	43400000	chr11	100000
chr11	43400001	48800000	Potocki-Shaffer	100000
chr11	48800001	114600000	chr11	100000
chr11	114600001	135086622	Jacobsen	100000
chr12	1	33275309	chr12	100000
chr13	1	114364328	chr13	100000
chr14	1	107043718	chr14	100000
chr15	1	20500000	chr15	100000
chr15	24500001	28193120	Angelman/Prader-Willi	100000
chr15	28193121	101991189	chr15	100000
chr16	1	17000000	chr16	100000
chr16	17000001	90338345	chr16	100000
chr17	1	3400000	Miller-Dieker	100000
chr17	3400001	16869758	chr17	100000
chr17	16869759	20318836	Smith-Magenis	100000
chr17	20318837	27400000	chr17	100000
chr17	27400001	33500000	NF1-microdeletion	100000
chr17	33500001	83257441	chr17	100000
chr18	1	80373285	chr18	100000
chr19	1	58617616	chr19	100000
chr20	1	64444167	chr20	100000
chr21	1	46709983	chr21	100000
chr22	1	17400000	chr22	100000
chr22	19022279	21098156	DiGeorge	100000
chr22	21098157	50818468	chr22	100000


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
devtools::install_github("seqinfo/BINDEL")
```
## Reference creation
BINDEL requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires the existence of `coordinates.bed` file, which defines the subregions in the genome to analyse.
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
`bindel::infer_normality("path/bam.bam", "reference_location.tsv")` outputs by default three scientific files:
1. `.png` illustrating high-risk probability per each region and reference sample set info for the same regions. The direction of the triangle hints if the findings are duplications or deletions. Each finding should be double-checked with the region bin figure.
2. `.png` containing normalised Z-score for each bin in each subregion. These Z-scores are the basis of the high-risk probability calculation. These figures also illustrate the reference group Z-scores.
3. `.tsv` summary file for each subregion.
