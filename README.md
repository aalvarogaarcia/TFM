# Unsupervised Clustering

This repository contains a Shiny user interface and an R Markdown report for unsupervised clustering for endophenotype discovery in genetic and proteomics data. 
The project was developed as part of a Master's thesis in Bioinformatics and biostatistics from UOC-UB in collaboration with the Genomics Division of the Institute of Technology and Renewable Energy (https://www.iter.es/areas/area-genomica/).

## Table of contents

* [Objective](#objective)
* [Methodolgy](#methodolgy)
* [Target Audience](#target-audience)
* [Disclaimer](#disclaimer)
* [Next Steps](#next-steps)

## Objective

The main objective of this project was to develop a user-friendly tool for exploring and analyzing genetic and protein data using unsupervised clustering techniques. 
The tool should allow users to:

* Perform univariate and bivariate analysis of genetic and protein data.
* Conduct multivariate analysis to identify endotypes using various clustering algorithms.
* Visualize the results of the clustering analysis using different visualization techniques.


## Methodology

The project was developed using R and Shiny. The Shiny user interface provides a user-friendly way to interact with the data and perform the analysis. 
The R Markdown report provides a detailed overview of the methodology and results.

### Key Features:

* Interactive Shiny user interface: The Shiny user interface allows users to easily upload data, select clustering algorithms and visualization techniques, and view the results of the analysis.
* Comprehensive R Markdown report: The R Markdown report provides a detailed overview of the methodology, including data preprocessing, clustering analysis, and visualization.
* Multiple R Scripts: In these scripts a whole series of methods and functions are developed with different objectives for the analysis and visualization of data.

## Target Audience:

 This tool is intended for researchers and students in the fields of bioinformatics, genetics, and genomics who are interested in analyzing genetic and protein data using unsupervised clustering techniques.

## Disclaimer:

 This tool is for research purposes only and should not be used for clinical diagnosis or treatment.

## Next Steps:

 Future work could include:

* Extending the tool to support additional data types.
* Developing more advanced visualization techniques.
* Integrating the tool with other bioinformatics tools.

<hr>
<!-- ------------------ SECTION 12 ------------------ -->

<a name="Update logs"></a>
## Update logs

> December 29, 2023. Several updates follows: figure 3 is updated; the bioinformatic workflow is enriched with an addition step of de novo assembly with Unycicler (SPAdes) in the case of multiple sequence alignment failure; BED files for Influenza A (H1N1 and H3N2) and B have been updated (strands and coordinates for the oligos were updated from a previous version of these files); useful BED files for mosdepth are provided to compute mean coverage-per-region in each strain.

> September 29, 2023. This repository became fully public. Enjoy the reading! ;=)

> September 26, 2023. Updated many sections: bioinformatic pipeline, primer-schemes (required BED files for the pipelines), deposited sequences, Influenza virus A and B reference sequences, and other useful repositories with resources to study Influenza.
 
> July 26, 2023. Created the private version of this repository.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
