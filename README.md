# BinDel: software tool for detecting clinically significant microdeletions in low-coverage WGS-based NIPT samples

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

## Alignment/Mapping
BinDel requires `.bam` **GRCh38** alignment files which are **sorted** and **duplicate marked**.

## Usage
### Reference creation
BinDel requires the creation of a reference set file. A reference set file is a file that contains known euploid NIPT samples. The read counts of these samples are used to compare the sample of interest with the healthy reference group. The creation of the reference file requires known euploid NIPT samples in `.bam` format and
a file defining genomic bins to use in the analysis.

```R
# In R:
bindel::create_reference("path/folder/bams", "path/bins.bed", "name_of_the_output_reference_file")
```

<details><summary> Click here to see how to create a file that defines genomic bins.</summary>
<p>

Given a file `locations.info.tsv` (column `length` sets bin size, e.g in here bin size is 300k):
```TSV
chr	start	end	focus	length
chr1	1	23599999	chr1	300000
chr1	23600000	27600000	1p36	300000
chr1	27600001	248956422	chr1	300000
chr2	1	242193529	chr2	300000
chr3	1	192599999	chr3	300000
chr3	192600000	198295559	3q29	300000
chr4	1	4500000	Wolf-Hirschhorn	300000
chr4	4500001	190214555	chr4	300000
chr5	1	10000	chr5	300000
chr5	10001	12533192	Cri-du-chat	300000
chr5	12533193	181538259	chr5	300000
chr6	1	170805979	chr6	300000
chr7	1	72699999	chr7	300000
chr7	72700000	77900000	Williams-Beuren	300000
chr7	77900001	159345973	chr7	300000
chr8	1	116699999	chr8	300000
chr8	116700000	126300000	Langer-Giedion	300000
chr8	126300001	145138636	chr8	300000
chr9	1	138394717	chr9	300000
chr10	1	133797422	chr10	300000
chr11	1	114599999	chr11	300000
chr11	114600000	135086622	Jacobsen	300000
chr12	1	133275309	chr12	300000
chr13	1	114364328	chr13	300000
chr14	1	107043718	chr14	300000
chr15	1	22677344	chr15	300000
chr15	22677345	28193120	Angelman/Prader-Willi	300000
chr15	28193121	101991189	chr15	300000
chr16	1	90338345	chr16	300000
chr17	1	83257441	chr17	300000
chr18	1	80373285	chr18	300000
chr19	1	58617616	chr19	300000
chr20	1	80105	chr20	300000
chr20	80106	1311812	20p13del	300000
chr20	1311813	64444167	chr20	300000
chr21	1	46709983	chr21	300000
chr22	1	17400000	chr22	300000
chr22	19022279	21098156	DiGeorge	300000
chr22	21098157	50818468	chr22	300000

```
Run the following Python script:
```
# In Python 3.x:
python dividebins.py --infile locations.info.tsv --outfile bins.bed
```
The script creates the file `bins.bed`, which can be used in the reference file creation.

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
bindel::infer_normality("path/bam.bam", "path/reference.tsv")
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
