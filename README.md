# Deciphering Multi-way Interactions in the Human Genome

# Overview
The 4DNvestigator is a MATLAB toolbox that analyzes time-series genome-wide chromosome
conformation capture (Hi-C) and gene expression (RNA-seq) data.

Paper: in preparation

Availability: https://github.com/lindsly/4DNvestigator

Data Availability: [4DNvestigator Data](https://drive.google.com/drive/folders/1xVjX7yqiOIPV_IfVKVDJJ79Ee0xMHGr8?usp=sharing)

# Installation
Install the 4DNvestigator through one of the following methods. 
Both methods will create a directory named 4DNvestigator which will contain all the necessary files for execution.
- Download the toolbox as a .zip folder and extract the contents into your MATLAB directory of choice. 
- Clone this Git repository into your MATLAB directory of choice using MATLAB's built-in [Source Control](https://www.mathworks.com/help/matlab/matlab_prog/retrieve-from-git-repository.html)

After the 4DNvestigator directory (and its sub-directories) are added, all 
[example data](https://drive.google.com/drive/folders/1xVjX7yqiOIPV_IfVKVDJJ79Ee0xMHGr8?usp=sharing) 
must be downloaded to the folder "data\exampleData" in order to run all of the functions contained in the example script "ExampleScript.m"

# Hi-C and RNA-seq data types
The 4DNvestigator accepts the following input file formats:

|**Data Type**|**File Type**|**Program**|
|----|----|----|
|Hi-C|.hic|Juicer|
|RNA-seq|.genes.results|RSEM|

# Core Functions
- 4DN Feature Analyzer: This measures the amount of change in both genome
structure and function for specified genomic regions by mapping all time
points to a consistent low dimensional embedding, and quantifying the variance
of each loci within this space over time. Method specifics can be found in:
["Genome Architecture Mediates Transcriptional Control of Human Myogenic Reprogramming"](https://www.cell.com/iscience/fulltext/S2589-0042(18)30114-7)
- von Neuman Entropy: Measures the entropy ("uncertainty") of a
multivariate system. Uncertainty is related to stemness. Here, we use this
measure to determine stemness of Hi-C samples.
- Larntz-Perlman: Method for testing the equality of correlation
matrices. This is applied to Hi-C correlation matrices to determine the
significance of matrix differences.
- Network Analysis: Measures the overlapping degree and multiplex participation 
coefficient to evaluate the heterogeneity of nodal degrees in dynamic biological networks.

# Other Functions
- Chromatin partitioning: Partitioning of the genome into two distinct
groups based on either the Fiedler vector or principal component 1. This
partitioning corresponds to euchromatin and heterochromatin, or A/B
compartments.
- Differential expression: Differential expression measures the
significance of RNA-seq expression differences between samples.
- A/B switching: This function determines which genomic regions change
their chromatin structure from compartment "A" to compartment "B"

# Examples
- See "ExampleScript.m" to run each of the core functions of the 4DNvestigator 
with default parameters using provided example data.
- [GettingStarted_4DNvestigator](https://github.com/lindsly/4DNvestigator/blob/master/MATLAB_Documentation/GettingStarted_4DNvestigator.pdf)
provides further details for each core function and their respective output.

# Dependencies
- The 4DNvestigator requires that [java](https://www.java.com/en/download/help/download_options.xml) and [python](https://www.python.org/downloads/) are installed to extract data from .hic files.
  - We recommend that the ["requests"](https://realpython.com/python-requests/) library is installed for python, and that it is added to the MATLAB path with "setenv"
- The 4DNvestigator uses [juicertools](https://github.com/aidenlab/juicer) to read data from .hic files
- The 4DNvestigator requires the following MATLAB Toolboxes for proper functionality: [Statistics and Machine Learning](https://www.mathworks.com/products/statistics.html), [Computer Vision](https://www.mathworks.com/products/computer-vision.html), [Bioinformatics](https://www.mathworks.com/products/bioinfo.html), and [Image Processing](https://www.mathworks.com/products/image.html)

# System Requirements
We recommend the following system properties for all 4DNvestigator functionalities to run properly:
- Windows 10
- At least 12 Gb RAM 
  - This is highly dependent on the size of the data being analyzed (data used in "ExampleScript.m" peaks at 9.5 Gb RAM usage) 
- A current version of MATLAB (2019 or later)

*4DNvestigator was written and tested on a Windows 10 machine with an Intel Core i7-8700 CPU and 32 Gb RAM using MATLAB R2019b.*
