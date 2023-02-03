### Minimum Working Example of Global River Slope Processing for Review Process. ###
Daniel Scherer, 2023-02-03  
[![DOI](https://zenodo.org/badge/595097911.svg)](https://zenodo.org/badge/latestdoi/595097911)

Note, that the actual IRIS software depends on the DGFI-TUM internal Multi-Version-Altimetry (MVA) database to access the ICESat-2 data efficiently.
Due to its size, the MVA database can not be shared.  
Therefore, this repository includes an executable jupyter notebook (*IRIS_example.ipynb*) that computes the water surface slope for one exemplary reach for which the input data has already been extracted from MVA.

The most relevant parts of IRIS are lines 977 to 1518 in dgfi_if/SWORD_IF.py:
- *get_icesat2_data* gets the required ICESat-2 data for the current reach from the MVA database.
- *icesat2_pass_helper* does the outlier rejection and along-track slope estimation.
- *get_icesat2_slope* and *icesat2_slope_helper* do the across-track and combined slope estimation.

*main.py* contains an example script to compute the slope globally by iterating over all SWORD reaches.

*export.py* contains the scripts to export the results to the published file format.

We describe the processing detailly in https://doi.org/10.1029/2022WR032842.  
Please contact daniel.scherer@tum.de for any further questions.
