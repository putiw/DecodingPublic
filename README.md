<a name="readme-top"></a>

<br />
<div align="center">

<h2 align="center">Identifying cortical areas that underlie the transformation from 2D retinal to 3D head-centric motion signals</h2>

  <p align="center">
  
  [![Screen][screenshot]](https://github.com/putiw/DecodingPublic/blob/main/helper_functions/chart.png)
  
  
  </p>
</div>

<!-- TABLE OF CONTENTS -->

  <summary><strong>Table of Contents</strong></summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#download-dataset-from-openneuro">Download dataset from OpenNeuro</a></li>
        <li><a href="#setup-github-repositories">Setup GitHub repositories</a></li>
      </ul>
    </li>     
    <li>
      <a href="#pipelines">Pipelines</a>
      <ul>
        <li><a href="#dcm2bids">dcm2bids</a></li>
        <li><a href="#fmriprep">fMRIPrep</a></li>
        <li><a href="#extract-rois">Extract ROIs</a></li>
      </ul>
    </li>
    <li>
      <a href="#analysis">Analysis</a>
      <ul>
        <li><a href="#tafkap">TAFKAP</a></li>
        <li><a href="#searchlight">Searchlight</a></li>
      </ul>
    </li>
    <li>
      <a href="#visualization">Visualization</a>
    </li>
  </ol>


<br />       


<!-- GETTING STARTED -->
## Getting Started

<br />  


### Download dataset from OpenNeuro

<br />  

Dataset link: [OpenNeuro](https://openneuro.org/datasets/ds004443/download)

<br />   

Download with *nodejs*, *S3*, or *datalad*
```sh
openneuro download --draft ds004443 ds004443-download/
aws s3 sync --no-sign-request s3://openneuro.org/ds004443 ds004443-download/
datalad install https://github.com/OpenNeuroDatasets/ds004443.git
```

<br />   

Dataset is in BIDS format such that the folder **rawdata** contains the orginal raw NIfTI files for both anatomical and funcational files, and the folder **derivatives/fmriprep** contains files that have been processed through fMRIPrep (in T1w space and in fsaverage6 space). Read more about the output space [here](https://fmriprep.org/en/stable/spaces.html).
 
<br />  
 
### Setup GitHub repositories

<br />  

1. Go to Github directory
  ```sh
  cd ~/Documents/GitHub 
  ```
2. Clone the repositories
  ```sh
  git clone https://github.com/putiw/DecodingPublic.git
  git clone https://github.com/Rokers/TAFKAP.git
  git clone https://github.com/cvnlab/knkutils.git
  git clone https://github.com/cvnlab/cvncode.git
  git clone https://github.com/gllmflndn/gifti.git
  ``` 

<p align="right">(<a href="#readme-top">back to top</a>)</p>
<br />       
<!-- Pipline -->


<br />   


## Pipelines

<br />   


### *dcm2bids*

<br />   


&nbsp;&nbsp;&nbsp;&nbsp;In this step we converted data from raw DICOM images into NIfTI files and organized data in BIDS format using [*Dcm2Bids*](https://unfmontreal.github.io/Dcm2Bids/). 

&nbsp;&nbsp;&nbsp;&nbsp;The output of this step is provided on OpenNeuro and is used in the proceeding steps. (The DICOM images will be provided upon request).

<br />   

### *fMRIPrep*

<br />   

&nbsp;&nbsp;&nbsp;&nbsp;We processed raw NIfTI files (those within **rawdata** folder) using [*fMRIPrep*](https://fmriprep.org/en/stable/installation.html).

&nbsp;&nbsp;&nbsp;&nbsp;The output of this step can be found in the **derivatives/fmriprep** folder. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:memo: **Note:** specifically, we ran fmriprep with this additioanl flag: *--output-space T1w:res-native fsaverage6:den-41k*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;:memo: **Note:** this step also requires [freeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/rel7downloads).

<br />   

### Extract ROIs

<br />   


&nbsp;&nbsp;&nbsp;&nbsp;We extracted the region of interests (ROIs) using the [neuropythy](https://github.com/noahbenson/neuropythy) pipeline. 

&nbsp;&nbsp;&nbsp;&nbsp;* Example way of running neuropythy with [Docker](https://docs.docker.com/engine/install/).

```sh
docker run -ti --rm -v <path-to-SUBJECTS_DIR>:/subjects nben/neuropythy atlas --verbose <subject-ID> --volume-export
```

&nbsp;&nbsp;&nbsp;&nbsp;This step generates the ROI atlas file named ***wang15_mplbl.mgz*** in **derivatives/freesurfer/sub-xxxx/mri** folder.

&nbsp;&nbsp;&nbsp;&nbsp;:memo: **Note:** You can find the output of this step (i.e. pre-generated the ROIs files in appropriate resolution) in **derivatives/fmriprep/sub-xxxx/ses-01/anat/rois** folder. 

&nbsp;&nbsp;&nbsp;&nbsp;:memo: **Note:** Alternatively, you can run [extractRois.m](https://github.com/putiw/DecodingPublic/blob/master/extractRois.m) in the current repo to generate the same files.

&nbsp;&nbsp;&nbsp;&nbsp;:memo: **Note:** A better way: [atlasmgz](https://github.com/WinawerLab/atlasmgz) (the paper did not use this way)

<br />   

## Analysis

<br />   

### TAFKAP

<br />   

This analysis resulted in figure 3 to 6. 
See orginal paper [here](https://www.biorxiv.org/content/10.1101/2021.03.04.433946v1).
See the code for this analysis here: [run_TAFKAP](https://github.com/putiw/DecodingPublic/blob/master/run_TAFKAP.m).

The input to this code is the fmriprep processed NIfTI files in **derivatives/fmriprep/sub-xxxx/ses-xx/func** folder. 
The file name should contain "*space-T1w_desc-preproc_bold.nii.gz".

The only variable that needs user input is the *bidsDir*, which should be set to whichever the directory of the project is in.

The code will print out the decoding accuracy for each ROI and for each subject as the output in the command window. These stats will be enough for you to regenerate the analysis in the paper, if you wish, you can save the result in a folder in **derivatives** for visualization purpose later. 

You can play with different ways of denoise, detrend, and estimate response amplitude in the sub-function [load_vol.m](https://github.com/putiw/DecodingPublic/blob/master/helper_functions/load_vol.m)

<br />   

### Searchlight


<br />   

This analysis resulted in figure 7 and 8. 
See orginal paper [here](https://www.biorxiv.org/content/10.1101/2021.03.04.433946v1).
See the code for this analysis here: [run_Searchlight](https://github.com/putiw/DecodingPublic/blob/master/run_Searchlight.m).

The required input files to the searchlight analysis are the functional scans in fsaverage6 space in the **derivatives/fmriprep/sub-xxxx/func/** folder.
The file name should contain "*space-fsaverage6_hemi-R_bold.func.gii".

You will also need the inflate surface files in the **derivatives/freesurfer/fsaverage6/mri** folder.

All of the above can be found on OpenNeuro.

<br />  

The output is the 81924 by 1 matrix containing one decoding accuracy value for each vertex in the fsaverage6 space. This result can be then used in further analysis. 


You will need to manually define three directories. 
```
bidsDir = '/Volumes/Vision/MRI/DecodingPublic'; % project directory containing the folder **derivatives**.
fsDir = '/Applications/freesurfer/7.2.0'; % freesurfer directory, most likely in application folder.
gitDir = '~/Documents/GitHub'; % GitHub directory, where you installed all the repositories. 
set_up(bidsDir,gitDir,fsDir)
```

You can play with different ways of denoise, detrend, and estimate response amplitude in the sub-function [load_surf.m](https://github.com/putiw/DecodingPublic/blob/master/helper_functions/load_surf.m)

<br />  

## Visualization

See [plot_figures.m](https://github.com/putiw/DecodingPublic/blob/master/helper_functions/plot_figures.m)


<!-- CONTACT -->
## Contact

Puti Wen - pw1246@nyu.edu



<p align="right">(<a href="#readme-top">back to top</a>)</p>






<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[screenshot]: helper_functions/chart.png
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
