# pySpade: Single cell Perturbations - Analysis of Differential gene Expression

## Overview
_________
`pySpade` is a user friendly tool to perform the whole transcriptome analysis of single cell perturbation dataset. With the direct output of Cellranger, `pySpade` utilizes hypergeomtric test to analyze the whole transcriptome differential expression and generates hits table csv file. User can use the table to do downstream processing like generating Manhattan plots (tutorial includes). Currently we support human genome.     

## Requirement
_________
* Python (3.7 +)
* Numpy (1.21 +)
* Pandas (1.3.5 +)
* Scipy (1.6.2 +)

## Installation
________
`pySpade` can be installed with `pip`

```shell
pip install pySpade
```

## Usage
________
```
$pySpade
usage: pySpade [-h]  ...

pySpade 
Version: 0.0.2

optional arguments:
  -h, --help  show this help message and exit

functions:
  
    process   process mapping output and reformat for downstream analysis.
    explevel  check the average expression level of query genes in single cell matrix
    fc        check the fold change of sgrna
    DEobs     perform differential expression analysis of observed cells
    DErand    perform differential expression analysis of random selection background
    local     perform local hit analysis with observation data and random background
    global    perform global hit analysis with observation data and random background
    
```

* `process` : Process transcriptome output and sgrna output to remove experimental doublets sgrna outlier cells. Transcriptome matrix is from Cellranger output (outs folder). sgrna matrix column: cell barcodes consistent with transcriptome matrix, rows: sgrna sequence. The final output format is h5 file. The final output can be compressed to save disk space, but it may take more time to write the final output file.

* `explevel` : Check the average expression level of query genes in single cell matrix. Input: processed transcriptome matrix from the `process` output, query genes list has to be txt file, genes are seperated with new line.

* `fc` : Check the fold change of perturbed region and individual sgRNA for query region and gene. Input: processed transcriptome and sgrna matrix from the `process` output and sgrna dict file (perturbation region hg38 coordinates and the sgrna sequence targeting that region. Region and sgrnas separated by tab, and sgrnas separated by comma)

* `DEobs` : Perfrom the genome wide differential expression analysis of all the perturbation regions. Input: processed transcriptome and sgrna matrix from the `process` output and sgrna dict file (perturbation region hg38 coordinates and the sgrna sequence targeting that region. Region and sgrnas separated by tab, and  sgrnas separated by comma). Output: up regulation p-value, downregulation p-value, fold change(compare with all the otehr cells) and average cpm.

* `DErand` : Perfrom the genome wide differential expression analysis of 1000 random selection cells. There are two options for random selection: all cells with equal probability or probability based on sgrna number in the cells. User should specify the cell number to select randomly. It is recommended with either exact cell number or bins (with large amount of perturbation experiment in order to reduce computational overhead).

* `local` : Using the observation p value and randomization bavckground p value to calculate the adjusted p value based on Gaussian distribution approximation. Local hits calculation includes the genes within plus and minus 2 Mb of the perturbation region. The output is a csv file with all hits information.

* `global` : Using the observation p value and randomization background p value to calculate the adjusted p value based on Gaussian distribution approximation. The output is a csv file with all hits information. Cutoff for p value is -1 (ln), fold change is 10%. Users can apply more stringent cutoff for their own purpose. 

## Contributors
_______
* First Author: Yihan Wang `Yihan.Wang@UTSouthwestern.edu`
* Corresponding Author: Gary Hon `Gary.Hon@UTSouthwestern.edu`