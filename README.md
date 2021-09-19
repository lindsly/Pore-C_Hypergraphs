# Deciphering Multi-way Interactions in the Human Genome

# Overview
All code to perform analysis in the manuscript, Deciphering Multi-way Interactions in the Human Genome, can be found here. The code is organized into separate folders which correspond to each figure in the main manuscript and supplemental materials.

Availability: https://github.com/lindsly/Pore-C_Hypergraphs

Data Availability: [PoreC Data](https://drive.google.com/drive/folders/1xVjX7yqiOIPV_IfVKVDJJ79Ee0xMHGr8?usp=sharing)

# Installation
Install the Pore-C_Hypergraphs code through one of the following methods. 
Both methods will create a directory named Pore-C_Hypergraphs which will contain all the necessary files for execution.
- Download the toolbox as a .zip folder and extract the contents into your MATLAB directory of choice. 
- Clone this Git repository into your MATLAB directory of choice using MATLAB's built-in [Source Control](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html).

After the Pore-C_Hypergraphs directory (and its sub-directories) are added, 
[data](https://drive.google.com/drive/folders/1xVjX7yqiOIPV_IfVKVDJJ79Ee0xMHGr8?usp=sharing) 
must be downloaded to the folder that is responsible for creating each figure.


# Code Organization
- Figure 2: Incidence matrix construction of a region in Chromosome 22 from fibroblasts. Contact frequency matrices were constructed by separating all multi-way contacts within this region of Chromosome 22 into their pairwise combinations. TADs are computed from the pairwise contacts. Multi-way contacts in this figure were determined in 100 kb resolution after noise reduction, originally derived from read-level multi-way contacts. (This code also applies to Supplemental Figure 1 for B lymphocytes).
- Figure 3: Incidence  matrix  construction  of  Chromosome  22  in fibroblasts. Frequencies of Pore-C contacts in Chromosome 22 can be shown in a bar plot according to the order of contact. Incidence matrix  construction  of  the  inter-chromosomal multi-way  contacts  between  Chromosome  20  and  Chromosome  22 in 1 Mb resolution.  Within this figure, all data are from one fibroblast sequencing run and multi-way contacts were determined after noise reduction at 1 Mb or 100 kb resolution accordingly.
- Figure 4: Incidence matrix construction of the top 10 most common multi-way contacts per chromosome. Matrices are constructed at 25 Mb resolution for both fibroblasts and B lymphocytes. Specifically, 5 intra-chromosomal and 5 inter-chromosomal multi-way contacts are identified for each chromosome with no repeated contacts. If 5 unique intra-chromosomal multi-way contacts are not possible in a chromosome, they are supplemented with additional inter-chromosomal contacts. Multi-way contacts were determined in 25 Mb resolution after noise reduction.
- Figure 5: The most common 2-way, 3-way, 4-way, and 5-way inter-chromosome combinations for each chromosome are calculated for fibroblasts and B lymphocytes. Inter-chromosomal combinations are determined using 25 Mb resolution multi-way contacts after noise reduction and are normalized by chromosome length. Here we only consider unique chromosome instances (i.e. multiple loci in a single chromosome are ignored).
- Figure 6: A 5 kb region before and after each locus in a Pore-C read is queried for chromatin accessibility and RNA Pol II binding (ATAC-seq and ChIP-seq, respectively). Multi-way contacts between accessible loci that have at least 1 instance of Pol II binding are indicative of potential transcription clusters. Gene expression and transcription factor binding sites are integrated to determine potential coexpression and coregulation within multi-way contacts with multiple genes. Transcription factor binding sites are queried +/- 5 kb from the gene’s transcription start site.
- Figure 7: Examples of potential transcription clusters are determined for fibroblasts and B lymphocytes. Multi-way contacts used for fibroblasts include all experiments. Examples were selected from the set of multi-way contacts summarized derived in the code for Figure 6.
- Figure S2: Calculates entropy of intra-chromosomal genomic hypergraph for fibroblast and B lymphocytes.
- Figure S3: Calculates hypergraph distance between two genome-wide hypergraphs derived from fibroblast and B lymphocytes. The background distribution is formed by measuring the hypergraph distances between the hypergraph derived from fibroblast and random hypergraphs. Also calculates hypergraph distances between intra-chromosomal genomic hypergraphs between fibroblasts and B lymphocytes. 

<!--# Dependencies
- The 4DNvestigator requires the following MATLAB Toolboxes for proper functionality: [Statistics and Machine Learning](https://www.mathworks.com/products/statistics.html), [Computer Vision](https://www.mathworks.com/products/computer-vision.html), [Bioinformatics](https://www.mathworks.com/products/bioinfo.html), and [Image Processing](https://www.mathworks.com/products/image.html)-->

# System Requirements
We recommend the following system properties for all code to run properly:
- Windows 10
- At least 16 Gb RAM 
  - This is highly dependent on the size of the data being analyzed
- A current version of MATLAB (2019 or later)

*This code was written and tested on a Windows 10 machine with an Intel Core i7-8700 CPU and 32 Gb RAM using MATLAB R2019b.*
