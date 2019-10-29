# DEM processing

### Filtering 
![Filtered DEM](images/filteredDEM.png)

### Cropping
![Cropped DEM](images/croppedDEM.png)

### Interpolating
![Interpolated DEM](images/interpolatedDEM.png)

### Repairing
![Repaired DEM](images/repairedDEM.png)

### Input parameters
| Parameters     | Value       | Description  |
| ------------- |-------------| ----------|
| Filtering   | *k1*, *k2* | Wavenumbers which define the low-pass filter |
| Cropping E   | *easting1*, *easting2*            | Only vertical or 3 components |
| Cropping N   | *northing1*, *northing2*            | Only vertical or 3 components |
| Interpolating | *dx*,*dy* | Site effect correction of signals |
| Repairing  | {*True*, *False*} | Length (seconds) of sliding time window |
| Input name  | string | Step size (seconds) of sliding time window |
| Output name | string | Corner values (Hertz) of frequency band which is used for localization  |

Generally it is best to filter before cropping in order to avoid filter artefacts in the domain of interest. However, for large DEM files it might be better to first crop the domain to reduce computational costs. 
