# BINDEL

The scientific software focuses on detecting rare (occurring in five or fewer people in 10,000) clinically relevant pathogenic microdeletions from low-coverage NIPT WGS data. 
However, the software is not limited to microdeletion detection and is developed with the idea of detecting any difference from the reference set. Detection possibility includes the detection of full chromosome aneuploidies or monosomies. When a female fetus reference set is used, then the software supports 45,X detection.

For examples and usage, please consult with the Cromwell workflow in the corresponding GitHub [repository](https://github.com/seqinfo/PPDxWorkflow).

# Manual

## Alignment
BINDEL requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Coordinates
BINDEL requires a `.bed` file with predefined coordinates. Coordinate file describes:
* Chromosomes (column `chr`)
* *bins* and bin lengths (bins can have varying length) which are defined by columns `start` and `end`.
* Names of the subregions of interest to analyse (column `focus`). Can be a whole chromosome or a subregion.

**Note:** Columns `chr`, `start` and `end` must uniquely define each region, e.g `.bed` file must not contain duplicates. Column `focus` is the name of the region of interest, which means that this column is used for grouping bins. **Having duplicates in .bed leads to anomalies in final high-risk probabilites**.

**Note 2:"** GC% correct depends on the number of regions of interest. E.g if only, for example chromosome 2 is in the analysis, it can affect the risk scoring compared to the having all chromosomes in the analysis.

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
chr	start	end	focus
chr1	1	27600000	1p36
chr1	27600001	248956422	chr1
chr2	1	242193529	chr2
chr3	1	198295559	chr3
chr4	1	4500000	Wolf-Hirschhorn
chr4	4500001	190214555	chr4
chr5	1	48800000	Cri-du-chat
chr5	48800001	181538259	chr5
chr6	1	170805979	chr6
chr7	1	72700000	chr7
chr7	72700001	77900000	Williams-Beuren
chr7	77900001	159345973	chr7
chr8	1	116700000	chr8
chr8	116700001	126300000	Langer-Giedion
chr8	126300001	145138636	chr8
chr9	1	138394717	chr9
chr10	1	133797422	chr10
chr11	1	114600000	chr11
chr11	114600001	135086622	Jacobsen
chr12	1	33275309	chr12
chr13	1	114364328	chr13
chr14	1	107043718	chr14
chr15	1	20500000	chr15
chr15	24500001	27800000	PWS/AS
chr15	27800001	101991189	chr15
chr16	1	17000000	chr16
chr16	17000001	90338345	chr16
chr17	1	83257441	chr17
chr18	1	80373285	chr18
chr19	1	58617616	chr19
chr20	1	64444167	chr20
chr21	1	46709983	chr21
chr22	1	17400000	chr22
chr22	19100001	21000000	DiGeorge
chr22	21700001	50818468	chr22

```
Run the following Python script:
```python
bin_width = 100000

with open("locations.info.tsv", encoding = "UTF-8", mode = "r") as f, open("all_regions.bed",  encoding = "UTF-8", mode = "w") as out:
    header = f.readline()
    out.write(header.strip() + "\n")
    for line in f:
        line = line.strip().split("\t")
        chromosome, start, end, focus = line[0], int(line[1]), int(line[2]), line[3]        
    
        while start + bin_width < end:
            out.write(f"{chromosome}\t{start}\t{(min(start + bin_width - 1, end))}\t{focus}\n")
            start = start + bin_width
```

## Installation
```R
# In R:
install.packages("devtools")
devtools::install_github("seqinfo/BINDEL")
```
## Reference creation
```R
# In R:
bindel::create_reference("path/folder/bams", "path/coordinate_bed.bed", "name_of_the_output")
```
## Usage
```R
# In R:
bindel::infer_normality("path/bam.bam", "reference_location.tsv")
```
## Output
`bindel::infer_normality("path/bam.bam", "reference_location.tsv")` outputs three scientific files:
1. `.png` illustrating high risk probability per each region and reference info for same regions.
2. `.png` containing normalised Z-scores per bins.
3. `.tsv` summary file for each subregion.
## Algorithm
The algorithm applies for each bin GC% correct, normalises bins by bin length and sample read count, applies PCA-based normalisation, calculates Z-scores based on the reference, applies Z-score normalisation based on the region mean read count, calculates Mahalanobis distance from the reference and converts them via Chi-Square distribution to high-risk probabilites.
