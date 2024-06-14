 <a name="TFM"></a> 

# Unsupervised Clustering

This repository contains a Shiny [user interface](https://github.com/aalvarogaarcia/TFM/blob/main/UI/user-interface.md) and an [R Markdown](https://github.com/aalvarogaarcia/TFM/tree/main/RMD/RMD.md) report for unsupervised clustering for endophenotype discovery in genetic and proteomics data. 
The project was developed as part of a Master's thesis in Bioinformatics and biostatistics from UOC-UB in collaboration with the Genomics Division of the Institute of Technology and Renewable Energy (https://www.iter.es/areas/area-genomica/).

## Table of contents

* [Objective](#objective)
* [Methodology](#methodology)
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
<!-- ------------------ SECTION 7 ------------------ -->

<a name="References"></a>
## References

<ol>
<li><a href="https://www.who.int/publications/i/item/9789240018440">Genomic sequencing of SARS-CoV-2. A guide to implementation for maximum impact on public health</a>, WHO, January 8, 2021.</li>
<li><a href="https://apps.who.int/iris/handle/10665/3"> Report “Global genomic surveillance strategy for pathogens with pandemic and epidemic potential, 2022-2032”</a>. Ginebra, WHO, 2022.</li>
<li>Gohl DM, Garbe J, Grady P, et al. A rapid, cost-effective tailed amplicon method for sequencing SARS-CoV-2. <i>BMC Genomics</i>. 2020;21(1):863. Published 2020 Dec 4. <a href="https://doi.org/10.1186/s12864-020-07283-6">doi:10.1186/s12864-020-07283-6</a>.</li>
<li>Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M. Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. <i>PLoS One</i>. 2020;15(9):e0239403. Published 2020 Sep 18. <a href="https://doi.org/10.1371/journal.pone.0239403">doi:10.1371/journal.pone.0239403</a>.</li> 
<li>Koskela von Sydow A, Lindqvist CM, Asghar N, et al. Comparison of SARS-CoV-2 whole genome sequencing using tiled amplicon enrichment and bait hybridization. <i>Sci Rep</i>. 2023;13(1):6461. Published 2023 Apr 20. <a href="https://doi.org/10.1038%2Fs41598-023-33168-1">doi:10.1038/s41598-023-33168-1</a>.</li>
<li>Zhou B, Wentworth DE. Influenza A virus molecular virology techniques. <i>Methods Mol Biol.</i> 2012;865:175-192. doi:<a href="https://link.springer.com/protocol/10.1007/978-1-61779-621-0_11">10.1007/978-1-61779-621-0_11</a>.</li>
<li>Zhou B, Lin X, Wang W, et al. Universal influenza B virus genomic amplification facilitates sequencing, diagnostics, and reverse genetics. <i>J Clin Microbiol</i>. 2014;52(5):1330-1337. doi:<a href="https://doi.org/10.1128/jcm.03265-13">10.1128/JCM.03265-13</a>.</li>
<li>Ying Lin, Jeffrey Koble, Priyanka Prashar, Anita Pottekat, Christina Middle, Scott Kuersten, Michael Oberholzer, Robert Brazas, Darcy Whitlock, Robert Schlaberg, Gary P. Schroth. A sequencing and subtyping protocol for Influenza A and B viruses using Illumina® COVIDSeq™ Assay Kit. Protocols.io. doi:<a href="https://dx.doi.org/10.17504/protocols.io.n2bvj8mrxgk5/v1">dx.doi.org/10.17504/protocols.io.n2bvj8mrxgk5/v1
</a></li>
</ol>
  
<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>



<hr>
<!-- ------------------ SECTION 8 ------------------ -->

<a name="Acknowledgements"></a>
## Acknowledgements

This study has been funded by Cabildo Insular de Tenerife (CGIEU0000219140 and "_Apuestas científicas del ITER para colaborar en la lucha contra la COVID-19_"); by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, epidemiological surveillance based on massive sequencing, Personalized Medicine and Biotechnology (OA17/008 and OA23/043); and by the agreement between Consejería de Educación, Universidades, Cultura y Deportes del Gobierno de Canarias y Cabildo Insular de Tenerife, 2022-2025 (AC0000014697).

This study is also an activity within the project Consolidation of WGS and RT-PCR activities for SARS-CoV-2 in Spain towards sustainable use and integration of enhanced infrastructure and capacities in the RELECOV network (101113109 - RELECOV 2.0) of the EU4Health Programme (EU4H) by the European Health and Digital Executive Agency (HaDEA), under the coordination of Instituto de Salud Carlos III (ISCIII).

We acknowledge the researchers and their institutions who released influenza sequences through NCBI GenBank, GISAID, and ENA that are being used in our studies. 

We also thank the authors, the laboratories that originated and submitted the genetic sequences and the metadata for sharing their work, as shown on Nextstrain, and:
<ul>
  <li>Hadfield <i>et al</i>, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics (2018).</li>
  <li>Sagulenko <i>et al</i>, TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017).</li>
</ul>

<!-- We would like to acknowledge the contributions of several researchers and laboratories who share their preliminary results through the [Virological](https://virological.org/) website. -->

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>



<hr>
<!-- ------------------ SECTION 11 ------------------ -->

<a name="How-to-cite"></a>
## How to cite this work

> This work has not been publised yet. See 'License and Attribution' section to cite this repository.

> To use the deposited sequences at GISAID, please, acknowledge this work as recommended by GISAID. Find the 'GISAID acknowledge tables' <a href="https://github.com/genomicsITER/influenza/tree/main/sequences/acknowledgements">here</a>.

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 9 ------------------ -->

<a name="License and Attribution"></a>
## License and Attribution

This repository and data exports are released under the CC BY 4.0 license. Please acknowledge the authors, the originating and submitting laboratories for the genetic sequences and metadata, and the open source software used in this work (third-party copyrights and licenses may apply).

Please cite this repository as: _"Influenza repository of the Reference Laboratory for Epidemiological Surveillance of Pathogens in the Canary Islands (accessed on YYYY-MM-DD)"_. And do not forget to <a href="#How-to-cite">cite the paper</a> (see the section "How to cite" below) when it becomes available. 

<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>




<hr>
<!-- ------------------ SECTION 12 ------------------ -->

<a name="Update logs"></a>
## Update logs

> June 18, 2024. My TFM github repository became public!

> May 19, 2024. Working hard into my GH repository.

> May 3, 2024. The developed UI is now available here.
 
> March 1, 2024. Created the private version of this repository.

<p align="right">
  <a href="#TFM" title="Up">
    <img src="https://github.com/aalvarogaarcia/TFM/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
