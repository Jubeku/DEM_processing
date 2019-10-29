# DEM processing

This routine allows to process Digital Elevation Models (DEMs).
DEMs can be filtered, cropped, interpolated and repaired.

### Input parameters

The following parameters are defined in the file `input.txt`:

| Function      | Variable     | Type (unit) | Description  |
| ------------- |---------------|-------| ----------|
| Filtering     | *k1*, *k2*    | float (1/m) | Wavenumbers defining the low-pass filter (k1=0 if no filtering is required) |
| Cropping E    | *easting1*, *easting2*   | float (m) | New easting coordinates |
| Cropping N    | *northing1*, *northing2* | float (m) | New northing coordinates |
| Interpolating | *dx*, *dy*    |float (m) | New resolution step size |
| Repairing     | *repair_bool* |{True, False} | Specify if reparation is required |
| Input name    | *fileID*      |string | File name of DEM |
| Output name   | *outID*       |string | File name to save the processed DEM  |


## Example processing
The processing functions are presented in the following using a DEM-file of Piton de la Fournaise volcano, La Réunion.

### Filtering
This functions allows to low-pass filter the DEM in the wavenumber domain. For this wavenumbers *k1* and *k2* are specified, defining the corner wavenumber and the maximum wavenumber of the filter, respectively. This means, that topographic variations of wavelengths shorter than 1/*k2* are filtered and $\lambda < $1/*k2* are completely removed.  
![Filtered DEM](images/filteredDEM.png)

### Cropping
![Cropped DEM](images/croppedDEM.png)

### Interpolating
![Interpolated DEM](images/interpolatedDEM.png)

### Repairing
![Repaired DEM](images/repairedDEM.png)

Generally it is best to filter before cropping in order to avoid filter artefacts in the domain of interest. However, for large DEM files it might be better to first crop the domain to reduce computational costs. 
