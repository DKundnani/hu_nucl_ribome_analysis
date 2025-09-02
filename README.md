
<h1 align="center">Human ribome Analysis - Storici & Jonoska Lab</h1>
This work amalgamates efforts of multiple computational biologists and mathematicians to dilineate riboncleotide (units of RNA) found to be embedded in human nuclear genomics DNA of multiple cell types. When ribonucleotides are embedded in the DNA, they are in the form of ribonucleoside mono phosphate or futher termed as rNMP(s).
<a name="readme-top"></a>
<!--

[![Commits][Commits-shield]][Commits-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Website][website-shield]][website-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
-->

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="##Installation">Installation</a></li>
      <ul>
        <li><a href="###Getting-the-code">Getting the code</a></li>
        <li><a href="###Creating-the-enviroment-with-required-dependencies">Required format of files</a></li>
      </ul>
    </li>
    <li><a href="##Usage">Breakdown of each folder as independent analysis</a></li>
    <li><a href="##Contact">Contact</a></li>
    <li><a href="##Citations">Citations</a></li>
  </ol>
</details>

<!-- Installation -->
## Installation
### Getting the code
Get development version from [GitHub](https://github.com/) with:
```sh
git clone https://github.com/DKundnani/hu_nucl_ribome_analysis
```

### Required format of files
* Input files (bed) containing single nucleotide locations, mainly for rNMP data. [link for GEOdataset] (another single nucleotide data can also be experimented on!)
* Reference genome files (.fa and .fai) of the organism being used(Also used to generate bed files) [Link for genome browser]
* Range files in preferably bed format. [Link for other location that has these files]

<!-- USAGE -->
## Breakdown of each folder as independent analysis
Please go into each repository to get usage and other details for each analysis.
* rNMP_EF: rNMP Enrichment Factor and Frequency calculation with rNMP EF heatmaps
* distribution_analysis: rNMP enrichment distribution around a fixed coordinate(s) in the genome
* rNMP randomized distrubtion and sequence association: Creating random rNMPs in the genome and histogram analysis, with getting mono and di-nucleotide correlations with rNMPs in discrete bins across the genome
* omics_processing: Standard pipelines used for RNA-seq, DNA-seq, bisulphite-seq and MACS2 peak calling in the study
* feature_correlations: Multi-omics associations for all rNMPs and different bases of rNMPs with expression and methylation data
* strand_bias: Study of strand bias near Transcription Start Site(s) and gene region
* agarose gel analysis: Distribution analysis of denatured DNA in presence and absence of RNase HII enzyme
* fragment size analysis: Simulation and estimation of rNMPs causing breaks in the agarose gel analysis


<!-- CONTACT -->
## Contact
Deepali L. Kundnani - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)    [![LinkedIn][linkedin-shield]][linkedin-url] 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Citations
Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!+
* <b> Human ribomes reveal DNA-embedded ribonucleotides as a new type of epigenetic mark. </b>
Deepali Lalchand Kundnani, Taehwan Yang, Tejasvi Channagiri, Penghao Xu, Yeunsoo Lee, Mo Sun, Francisco Martinez-Figueroa, Supreet Randhawa, Ashlesha Gogate, Youngkyu Jeon, Stefania Marsili, Gary Newnam, Yilin Lu, Vivian Park, Sijia Tao, Justin Ling, Raymond Schinazi, Zachary Pursell, Abdulmelik Mohammed, Patricia Opresko, Bret Freudenthal, Baek Kim, Soojin Yi, Nataša Jonoska, Francesca Storici, <i>  bioRxiv </i> 2025, [https://doi.org/10.1101/2025.06.27.661996](https://doi.org/10.1101/2025.06.27.661996)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[contributors-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[forks-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/forks
[stars-shield]: https://img.shields.io/github/stars/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[stars-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/stargazers
[issues-shield]: https://img.shields.io/github/issues/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[issues-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/issues
[license-shield]: https://img.shields.io/github/license/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[license-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/deepalik
[product-screenshot]: images/screenshot.png
[commits-url]: https://github.com/DKundnani/hu_nucl_ribome_analysis/pulse
[commits-shield]: https://img.shields.io/github/commit-activity/t/DKundnani/hu_nucl_ribome_analysis?style=for-the-badge
[website-shield]: https://img.shields.io/website?url=http%3A%2F%2Fdkundnani.bio%2F&style=for-the-badge
[website-url]:http://dkundnani.bio/ 
