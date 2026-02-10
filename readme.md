## Code for my paper on the Sextans A PAH detection with JWST photometry

ADS link: https://ui.adsabs.harvard.edu/abs/2025arXiv251204060T/abstract

These scripts are sadly not well documented, but may help you get started. Descriptions of each folder:

- data reduction: contains the pipeline processing of the JWST imaging data, including NIRCam amd MIRI
- psf-match: relatively stand alone code to produce convolution kernels for JWST images. Implements the Aniano+ 2011 method for kernel generation
- models: applies synthetic photometry to Draine+ 2021 models and PDRS4All spectra. Identified the "k" values used in continuum subtraction. 
- con sub: bulk of the analysis code and generates most of the plots in the paper. 

 
