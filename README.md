# DEM processing

### Input parameters

The parameters are defined in the file `input.txt`.

| Function      | Variable     | Type (unit) | Description  |
| ------------- |---------------|-------| ----------|
| Filtering     | *k1*, *k2*    | float (1/m) | Wavenumbers defining the low-pass filter |
| Cropping E    | *easting1*, *easting2*   | float (m) | New easting coordinates |
| Cropping N    | *northing1*, *northing2* | float (m) | New northing coordinates |
| Interpolating | *dx*, *dy*    |float (m) | New resolution step size |
| Repairing     | *repair_bool* |{True, False} | Specify if reparation is wished |
| Input name    | *fileID*      |string | File name of DEM |
| Output name   | *outID*       |string | File name to save the processed DEM  |

### Filtering 
![Filtered DEM](images/filteredDEM.png)

### Cropping
![Cropped DEM](images/croppedDEM.png)

### Interpolating
![Interpolated DEM](images/interpolatedDEM.png)

### Repairing
![Repaired DEM](images/repairedDEM.png)

Generally it is best to filter before cropping in order to avoid filter artefacts in the domain of interest. However, for large DEM files it might be better to first crop the domain to reduce computational costs. 
