
<h1 align="center">RPA-wrapper</h1>
<a href="https://github.com/xph9876/RibosePreferenceAnalysis">Ribose Preferred Analysis</a> Wrapper Scripts

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>


[![Commits][Commits-shield]][Commits-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Website][website-shield]][website-url]
[![Issues][issues-shield]][issues-url]
[![License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="##Installation">Installation</a></li>
      <ul>
        <li><a href="###Getting-the-code">Getting the code</a></li>
        <li><a href="###Creating-the-enviroment-with-required-dependencies">Creating the enviroment with required dependencies</a></li>
        <li><a href="###Additional-Dependencies">Additional Dependencies</a></li>
      </ul>
    </li>
    <li><a href="##Usage">Usage</a></li>
      <ul>
        <li><a href="###Defining-variables">Defining variables</a></li>
        <li><a href="Initializing-functions-and-activating-enrviroment">Initializing functions and activating enrviroment</a></li>
        <li><a href="###Running-the-Code">Running the Code</a></li>
       <li><a href="###Statistical-test-for-preference-and-genotype-comparisons">Statistical test for preference and genotype comparisons</a></li>
        <li><a href="###Generating-Stacked-barplots-for-hotspots">Generating Stacked barplots for hotspots</a></li>
      </ul>
    <li><a href="##Contributing">Contributing</a></li>
    <li><a href="##License">License</a></li>
    <li><a href="##Contact">Contact</a></li>
    <li><a href="##Citations">Citations</a></li>
  </ol>
</details>

<!-- Installation -->
## Installation
### Getting the code
The development version from [GitHub](https://github.com/) with:
```sh
git clone https://github.com/xph9876/RibosePreferenceAnalysis.git
git clone https://github.com/DKundnani/RPA-wrapper.git
```
### Creating the enviroment with required dependencies
```sh
conda env create --name RPAwrapper_env --file /RPA-wrapper/env.yml
```
### Additional Dependencies
* Input files (bed) containing single nucleotide locations, mainly for rNMP data. (another single nucleotide data can also be experimented on!)
* Reference genome files (.fa and .fai) of the organism being used(Also used to generate bed files)
* ranges/bed file fo the genome locations to be analyzed and for which background frequency will be calculated as well.
* order file, example in repository

<!-- USAGE -->
## Usage
### Defining variables
```bash
scripts='path/to/RibosePreferenceAnalysis/' #location of RPA repository
ref='path/to/reference/sacCer2/sacCer2.fa' #Reference Fast file
range='path/to/ranges/sacCer2/nucl.bed' #Example ranges of nuclear genome of sacCer2
#range='path/to/chrM.bed' 
bed='path/to/bed/' #Folder of bed files
order='path/to/order' #sample file share in the RPA-wrapper repository
```
### Initializing functions and activating enrviroment
```bash
conda activate RPAwrapper_env #activating enviroment
source path/to/RPA-wrapper/Heatmapwrapper.sh
```
### Running the Code for both strands
```bash
mkdir heatmaps; cd heatmaps #make the output directory and run the code from it
#Generating Background frequency of the genome
bg_freq $scripts $ref $range #one time for each range
sample_freq $scripts $ref $range $bed #Generating frequency of libraries/samples
norm_freq $scripts $ref $range $bed #Normalizing sample frequency to genome frequency
resort_plot $scripts $ref $range $bed $order #resort the matrices as per order file and hence the heatmaps
rm sample_freq/$(basename $range .bed)/*bed #Cleanup
```

### Running the Code for single strand - gives out same and opposite for the range with 6th column having strand information
```bash
mkdir heatmaps; cd heatmaps #make the output directory and run the code from it
#Generating Background frequency of the genome
bg_freq_ss $scripts $ref $range ; bg_freq_os $scripts $ref $range
sample_freq_ss $scripts $ref $range $bed ; sample_freq_os $scripts $ref $range $bed
norm_freq_ss $scripts $ref $range $bed ; norm_freq_os $scripts $ref $range $bed
resort_plot_ss $scripts $ref $range $bed $order ; resort_plot_os $scripts $ref $range $bed $order
rm sample_freq/$(basename $range .bed)/*bed #Cleanup
```

### Statistical test for preference and genotype comparisons
```bash
mww $scripts $ref $range $bed $order
```

### Generating Stacked barplots for hotspots
```bash
Rscript path/to/RPA-wrapper/comp.R -m sorted_chrM_mono_0 #hotspot composition files usually contain on entry for every genotype.
```


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GNU GPL3 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact
Deepali L. Kundnani - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)    [![LinkedIn][linkedin-shield]][linkedin-url] 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->
## Citations
Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!+
* <b> Distinct features of ribonucleotides within genomic DNA in Aicardi-Goutières syndrome (AGS)-ortholog mutants of <i>Saccharomyces cerevisiae</i> </b>
Deepali L. Kundnani, Taehwan Yang, Alli L. Gombolay, Kuntal Mukherjee, Gary Newnam, Chance Meers, Zeel H. Mehta, Celine Mouawad, Francesca Storici
bioRxiv 2023.10.02.560505; doi:[https://doi.org/10.1101/2023.10.02.560505]( https://doi.org/10.1101/2023.10.02.560505)
* Kundnani, D. (2024). rNMP_hotspots:2.0.0 (2.0.0). Zenodo.  [https://doi.org/10.5281/zenodo.8152090](https://doi.org/10.5281/zenodo.8152090) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8152090.svg)](https://doi.org/10.5281/zenodo.8152090)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/DKundnani/RPA-wrapper?style=for-the-badge
[contributors-url]: https://github.com/DKundnani/RPA-wrapper/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/DKundnani/RPA-wrapper?style=for-the-badge
[forks-url]: https://github.com/DKundnani/RPA-wrapper/forks
[stars-shield]: https://img.shields.io/github/stars/DKundnani/RPA-wrapper?style=for-the-badge
[stars-url]: https://github.com/DKundnani/RPA-wrapper/stargazers
[issues-shield]: https://img.shields.io/github/issues/DKundnani/RPA-wrapper?style=for-the-badge
[issues-url]: https://github.com/DKundnani/RPA-wrapper/issues
[license-shield]: https://img.shields.io/github/license/DKundnani/RPA-wrapper?style=for-the-badge
[license-url]: https://github.com/DKundnani/RPA-wrapper/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/deepalik
[product-screenshot]: images/screenshot.png
[commits-url]: https://github.com/DKundnani/RPA-wrapper/pulse
[commits-shield]: https://img.shields.io/github/commit-activity/t/DKundnani/RPA-wrapper?style=for-the-badge
[website-shield]: https://img.shields.io/website?url=http%3A%2F%2Fdkundnani.bio%2F&style=for-the-badge
[website-url]:http://dkundnani.bio/ 
